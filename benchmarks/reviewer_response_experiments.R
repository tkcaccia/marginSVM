#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialGraphRefine)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(mclust)
})

out_dir <- file.path("benchmarks", "results", "reviewer_response")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
quick <- identical(Sys.getenv("SPATIAL_REFINE_QUICK"), "1")
reps <- if (quick) 3L else 10L
n_main <- if (quick) 6000L else 12000L
benchmark_workers <- max(1L, as.integer(Sys.getenv("SPATIAL_REFINE_BENCHMARK_WORKERS", "4")))

theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank()))

variants <- list(
  "Spatial graph" = list(),
  "C++ kNN vote" = list(weighted = FALSE, iterations = 1L, consensus = 0.5,
                          preserve = 0, margin = 0, current_support = 1),
  "Potts-like ICM" = list(weighted = FALSE, iterations = 8L, consensus = 0.5,
                           preserve = 0.35, margin = 0, current_support = 1)
)
weighted_one_pass <- list(weighted = TRUE, iterations = 1L, consensus = 0.5,
                          preserve = 0, margin = 0, current_support = 1)

safe_mean <- function(x) if (length(x)) mean(x) else NA_real_

calibration_error <- function(confidence, correct, bins = 10L) {
  if (is.null(confidence)) return(NA_real_)
  groups <- cut(confidence, breaks = seq(0, 1, length.out = bins + 1L), include.lowest = TRUE)
  tab <- split(seq_along(correct), groups, drop = TRUE)
  sum(vapply(tab, function(i) length(i) / length(correct) * abs(mean(confidence[i]) - mean(correct[i])), numeric(1)))
}

score <- function(pred, sim, elapsed, method) {
  truth <- sim$truth
  initial <- sim$labels
  recalls <- vapply(levels(truth), function(level) safe_mean(pred[truth == level] == level), numeric(1))
  smallest <- names(which.min(table(truth)))
  changed <- pred != initial
  initially_wrong <- initial != truth
  initially_correct <- !initially_wrong
  confidence <- attr(pred, "confidence")
  correct <- pred == truth
  prevalence_truth <- prop.table(table(truth))
  prevalence_pred <- prop.table(table(factor(pred, levels = levels(truth))))
  tibble(
    method = method,
    accuracy = mean(correct),
    ari = mclust::adjustedRandIndex(pred, truth),
    macro_recall = mean(recalls, na.rm = TRUE),
    worst_recall = min(recalls, na.rm = TRUE),
    smallest_recall = unname(recalls[smallest]),
    boundary_005 = safe_mean(correct[sim$boundary_proximity <= 0.05]),
    boundary_010 = safe_mean(correct[sim$boundary_proximity <= 0.10]),
    boundary_020 = safe_mean(correct[sim$boundary_proximity <= 0.20]),
    interior_accuracy = safe_mean(correct[sim$boundary_proximity > 0.20]),
    changed_precision = safe_mean(correct[changed]),
    correction_recall = safe_mean(correct[initially_wrong]),
    damage_rate = safe_mean(!correct[initially_correct]),
    changed_fraction = mean(changed),
    area_l1 = sum(abs(prevalence_truth - prevalence_pred)),
    ece = calibration_error(confidence, correct),
    brier = if (is.null(confidence)) NA_real_ else mean((confidence - correct)^2),
    seconds = unname(elapsed["elapsed"])
  )
}

run_methods <- function(sim, methods = variants) {
  initial <- score(sim$labels, sim, c(elapsed = 0), "Initial")
  rows <- lapply(names(methods), function(method) {
    timing <- system.time(pred <- refine_spatial_clusters(sim$xy, sim$labels, sim$samples,
                                                           control = methods[[method]]))
    score(pred, sim, timing, method)
  })
  timing <- system.time(spagcn <- SpatialGraphRefine:::.refine_published_labels(
    sim$xy, sim$labels, sim$samples, method = "spagcn", neighbors = 6L
  ))
  rows[[length(rows) + 1L]] <- score(spagcn, sim, timing, "SpaGCN refine")
  timing <- system.time(graphst <- SpatialGraphRefine:::.refine_published_labels(
    sim$xy, sim$labels, sim$samples, method = "graphst", neighbors = 50L
  ))
  rows[[length(rows) + 1L]] <- score(graphst, sim, timing, "GraphST refine")
  bind_rows(c(list(initial), rows))
}

message("Running frozen-default robustness experiment")
robust_design <- crossing(
  pattern = c("jagged_stripes", "wavy_layers", "rings", "spiral", "branching",
              "lobes", "islands", "disconnected", "thin_layers", "intermixed"),
  noise = if (quick) c(0.15, 0.30) else c(0.05, 0.15, 0.25, 0.40),
  noise_type = if (quick) c("random", "patch") else c("random", "boundary", "patch", "region"),
  replicate = seq_len(reps)
)

run_robust_scenario <- function(i) {
  d <- robust_design[i, ]
  sim <- simulate_spatial_domains(
    n_main, d$pattern, k = 5L, noise = d$noise, noise_type = d$noise_type,
    seed = 100000L + i
  )
  metrics <- run_methods(sim)
  bind_cols(d[rep(1L, nrow(metrics)), ], metrics)
}
robust_results <- if (.Platform$OS.type != "windows" && benchmark_workers > 1L) {
  parallel::mclapply(seq_len(nrow(robust_design)), run_robust_scenario,
                     mc.cores = benchmark_workers, mc.preschedule = TRUE)
} else {
  lapply(seq_len(nrow(robust_design)), run_robust_scenario)
}
robust_metrics <- bind_rows(robust_results)
write.csv(robust_metrics, file.path(out_dir, "robustness_metrics.csv"), row.names = FALSE)

message("Running ablations on held-out seeds")
ablations <- list(
  "Full method" = list(),
  "No distance weights" = list(weighted = FALSE),
  "No current-label prior" = list(preserve = 0),
  "No margin gate" = list(margin = 0),
  "No confidence gates" = list(consensus = 0.5, margin = 0),
  "No discordance gate" = list(current_support = 1),
  "One iteration" = list(iterations = 1L),
  "C++ kNN vote" = variants[["C++ kNN vote"]],
  "Weighted one-pass" = weighted_one_pass
)
ablation_design <- crossing(
  pattern = c("jagged_stripes", "islands", "thin_layers", "intermixed"),
  noise_type = c("random", "boundary", "patch"),
  replicate = seq_len(reps)
)
ablation_results <- vector("list", nrow(ablation_design))
for (i in seq_len(nrow(ablation_design))) {
  d <- ablation_design[i, ]
  sim <- simulate_spatial_domains(n_main, d$pattern, noise = 0.25, noise_type = d$noise_type,
                                  feature_scale = if (d$pattern %in% c("islands", "thin_layers")) 0.7 else 1,
                                  seed = 200000L + i)
  rows <- lapply(names(ablations), function(method) {
    timing <- system.time(pred <- refine_spatial_clusters(sim$xy, sim$labels, control = ablations[[method]]))
    score(pred, sim, timing, method)
  })
  ablation_results[[i]] <- bind_cols(d[rep(1L, length(rows)), ], bind_rows(rows))
}
ablation_metrics <- bind_rows(ablation_results)
write.csv(ablation_metrics, file.path(out_dir, "ablation_metrics.csv"), row.names = FALSE)

message("Running one-factor sensitivity analysis")
sensitivity_settings <- bind_rows(
  tibble(parameter = "neighbors", value = c(7, 11, 15, 21, 31, 45, 61),
         control = lapply(value, function(x) list(neighbors = as.integer(x)))),
  tibble(parameter = "consensus", value = c(0.50, 0.56, 0.65, 0.75),
         control = lapply(value, function(x) list(consensus = x))),
  tibble(parameter = "preserve", value = c(0, 0.06, 0.12, 0.25, 0.50),
         control = lapply(value, function(x) list(preserve = x))),
  tibble(parameter = "margin", value = c(0, 0.05, 0.10, 0.20),
         control = lapply(value, function(x) list(margin = x))),
  tibble(parameter = "iterations", value = c(1, 2, 3, 5),
         control = lapply(value, function(x) list(iterations = as.integer(x)))),
  tibble(parameter = "current_support", value = c(0.10, 0.15, 0.20, 0.25, 0.35, 1.00),
         control = lapply(value, function(x) list(current_support = x)))
)
sensitivity_results <- list()
counter <- 0L
for (pattern in c("jagged_stripes", "islands", "thin_layers", "intermixed")) {
  for (replicate in seq_len(if (quick) 2L else 5L)) {
    sim <- simulate_spatial_domains(n_main, pattern, noise = 0.25, noise_type = "patch",
                                    feature_scale = if (pattern %in% c("islands", "thin_layers")) 0.7 else 1,
                                    seed = 300000L + counter + replicate)
    for (i in seq_len(nrow(sensitivity_settings))) {
      counter <- counter + 1L
      timing <- system.time(pred <- refine_spatial_clusters(sim$xy, sim$labels,
                                                             control = sensitivity_settings$control[[i]]))
      sensitivity_results[[counter]] <- bind_cols(
        tibble(pattern = pattern, replicate = replicate,
               parameter = sensitivity_settings$parameter[i], value = sensitivity_settings$value[i]),
        score(pred, sim, timing, "Spatial graph")
      )
    }
  }
}
sensitivity_metrics <- bind_rows(sensitivity_results)
write.csv(sensitivity_metrics, file.path(out_dir, "sensitivity_metrics.csv"), row.names = FALSE)

message("Running rare and thin-domain survival experiment")
feature_design <- crossing(
  pattern = c("islands", "thin_layers"),
  feature_scale = c(0.35, 0.50, 0.70, 1.00, 1.40),
  noise = c(0.15, 0.25),
  replicate = seq_len(reps)
)
feature_results <- vector("list", nrow(feature_design))
for (i in seq_len(nrow(feature_design))) {
  d <- feature_design[i, ]
  sim <- simulate_spatial_domains(if (quick) 10000L else 20000L, d$pattern, noise = d$noise,
                                  noise_type = "random", feature_scale = d$feature_scale,
                                  seed = 400000L + i)
  metrics <- run_methods(sim)
  feature_results[[i]] <- bind_cols(d[rep(1L, nrow(metrics)), ], metrics)
}
feature_metrics <- bind_rows(feature_results)
write.csv(feature_metrics, file.path(out_dir, "feature_survival_metrics.csv"), row.names = FALSE)

message("Running domain-count experiment")
count_design <- crossing(
  pattern = c("jagged_stripes", "rings", "branching"),
  clusters = c(2L, 3L, 5L, 8L, 12L),
  replicate = seq_len(if (quick) 2L else 8L)
)
count_results <- vector("list", nrow(count_design))
for (i in seq_len(nrow(count_design))) {
  d <- count_design[i, ]
  sim <- simulate_spatial_domains(n_main, d$pattern, k = d$clusters, noise = 0.25,
                                  noise_type = "random", seed = 500000L + i)
  metrics <- run_methods(sim)
  count_results[[i]] <- bind_cols(d[rep(1L, nrow(metrics)), ], metrics)
}
count_metrics <- bind_rows(count_results)
write.csv(count_metrics, file.path(out_dir, "domain_count_metrics.csv"), row.names = FALSE)

message("Running clean-label negative controls")
negative_results <- list()
counter <- 0L
for (pattern in c("jagged_stripes", "thin_layers", "intermixed")) {
  for (replicate in seq_len(reps)) {
    counter <- counter + 1L
    sim <- simulate_spatial_domains(n_main, pattern, noise = 0, seed = 600000L + counter)
    negative_results[[counter]] <- bind_cols(
      tibble(pattern = pattern, replicate = replicate),
      run_methods(sim, variants[c("Spatial graph", "C++ kNN vote", "Potts-like ICM")])
    )
  }
}
negative_metrics <- bind_rows(negative_results)
write.csv(negative_metrics, file.path(out_dir, "negative_control_metrics.csv"), row.names = FALSE)

message("Running invariance and sample-stratification checks")
base <- simulate_spatial_domains(15000, "spiral", noise = 0.25, seed = 700001L)
base_pred <- refine_spatial_clusters(base$xy, base$labels)
angle <- 0.73
rotation <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2L)
transformed <- base$xy %*% rotation
transformed <- sweep(transformed, 2L, c(12.3, -7.8), "+")
rot_pred <- refine_spatial_clusters(transformed, base$labels)
permutation <- sample.int(nrow(base$xy))
permuted <- refine_spatial_clusters(base$xy[permutation, ], base$labels[permutation])
permuted_back <- permuted[base::order(permutation)]
label_map <- sample(levels(base$labels))
names(label_map) <- levels(base$labels)
relabeled <- factor(label_map[as.character(base$labels)], levels = label_map)
relabeled_pred <- refine_spatial_clusters(base$xy, relabeled)
restored <- factor(names(label_map)[match(as.character(relabeled_pred), label_map)], levels = levels(base$labels))

multi <- simulate_spatial_domains(30000, "wavy_layers", samples = 6, noise = 0.25, seed = 700002L)
multi_pred <- refine_spatial_clusters(multi$xy, multi$labels, multi$samples)
invariance <- tibble(
  test = c("rotation_translation", "row_order", "label_permutation", "multi_sample_accuracy"),
  agreement = c(mean(rot_pred == base_pred), mean(permuted_back == base_pred),
                mean(restored == base_pred), mean(multi_pred == multi$truth))
)
write.csv(invariance, file.path(out_dir, "invariance_metrics.csv"), row.names = FALSE)
write.csv(data.frame(sample = levels(multi$samples), neighbors = attr(multi_pred, "neighbors")),
          file.path(out_dir, "multisample_neighbors.csv"), row.names = FALSE)

message("Generating reviewer-response figures")
method_levels <- c("Initial", "Spatial graph", "SpaGCN refine", "GraphST refine",
                   "C++ kNN vote", "Potts-like ICM")
robust_summary <- robust_metrics |>
  mutate(method = factor(method, method_levels)) |>
  group_by(noise_type, noise, method) |>
  summarise(across(c(accuracy, boundary_010, smallest_recall, damage_rate),
                   function(x) mean(x, na.rm = TRUE)), .groups = "drop")

p_robust <- ggplot(robust_summary, aes(noise, accuracy, color = method)) +
  geom_line() + geom_point(size = 1.4) +
  facet_wrap(~noise_type) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Corrupted-label fraction", y = "Accuracy", color = "Method")
ggsave(file.path(out_dir, "robustness_by_error.png"), p_robust, width = 10, height = 5, dpi = 220)

p_boundary <- ggplot(robust_summary, aes(noise, boundary_010, color = method)) +
  geom_line() + geom_point(size = 1.4) +
  facet_wrap(~noise_type) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Corrupted-label fraction", y = "Accuracy near analytic boundary", color = "Method")
ggsave(file.path(out_dir, "analytic_boundary_accuracy.png"), p_boundary, width = 10, height = 5, dpi = 220)

ablation_summary <- ablation_metrics |>
  group_by(method) |>
  summarise(accuracy = mean(accuracy), boundary = mean(boundary_010),
            damage = mean(damage_rate), .groups = "drop") |>
  pivot_longer(c(accuracy, boundary, damage), names_to = "metric", values_to = "value")
p_ablation <- ggplot(ablation_summary, aes(reorder(method, value), value, fill = method)) +
  geom_col(show.legend = FALSE) + coord_flip() + facet_wrap(~metric, scales = "free_x") +
  labs(x = NULL, y = "Mean value")
ggsave(file.path(out_dir, "ablation_summary.png"), p_ablation, width = 9, height = 5, dpi = 220)

p_feature <- feature_metrics |>
  filter(method %in% c("Initial", "Spatial graph", "SpaGCN refine", "GraphST refine",
                       "C++ kNN vote", "Potts-like ICM")) |>
  group_by(pattern, feature_scale, noise, method) |>
  summarise(smallest_recall = mean(smallest_recall), .groups = "drop") |>
  ggplot(aes(feature_scale, smallest_recall, color = method)) +
  geom_line() + geom_point() + facet_grid(pattern ~ noise, labeller = label_both) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Relative feature width", y = "Smallest-domain recall", color = "Method")
ggsave(file.path(out_dir, "rare_thin_survival.png"), p_feature, width = 10, height = 6, dpi = 220)

p_sensitivity <- sensitivity_metrics |>
  group_by(pattern, parameter, value) |>
  summarise(accuracy = mean(accuracy), damage = mean(damage_rate), .groups = "drop") |>
  ggplot(aes(value, accuracy, color = pattern)) +
  geom_line() + geom_point(size = 1.2) + facet_wrap(~parameter, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) + labs(x = "Parameter value", y = "Accuracy", color = "Geometry")
ggsave(file.path(out_dir, "parameter_sensitivity.png"), p_sensitivity, width = 11, height = 6, dpi = 220)

calibration <- robust_metrics |>
  filter(method == "Spatial graph") |>
  summarise(ece = mean(ece, na.rm = TRUE), brier = mean(brier, na.rm = TRUE),
            changed_precision = mean(changed_precision, na.rm = TRUE),
            correction_recall = mean(correction_recall, na.rm = TRUE),
            damage_rate = mean(damage_rate, na.rm = TRUE))
write.csv(calibration, file.path(out_dir, "calibration_summary.csv"), row.names = FALSE)

message("Reviewer-response experiments complete: ", out_dir)
