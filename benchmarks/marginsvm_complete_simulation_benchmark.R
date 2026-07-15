#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialGraphRefine)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mclust)
})

out_dir <- file.path("benchmarks", "results", "marginsvm_complete_simulation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
quick <- identical(tolower(Sys.getenv("MARGINSVM_COMPLETE_QUICK", "false")), "true")
workers <- max(1L, as.integer(Sys.getenv("SPATIAL_REFINE_BENCHMARK_WORKERS", "4")))
n_geometry <- if (quick) 2500L else 6000L
replicates <- seq_len(if (quick) 1L else 3L)
patterns <- c("jagged_stripes", "wavy_layers", "rings", "spiral", "branching",
              "lobes", "islands", "disconnected", "thin_layers", "intermixed",
              "layers3d")
density_profiles <- c("uniform", "moderate", "strong", "extreme", "hotspot")
noise_types <- c("random", "boundary", "patch", "region")
method_levels <- c("Initial", "marginSVM", "marginSVM trust field", "Graph refinement",
                   "SpaGCN refine", "GraphST refine", "C++ kNN vote", "Potts-like ICM")

theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank()))

safe_mean <- function(x) if (length(x)) mean(x) else NA_real_

score_prediction <- function(pred, sim, elapsed, method) {
  correct <- pred == sim$truth
  initially_wrong <- sim$labels != sim$truth
  changed <- pred != sim$labels
  class_recall <- vapply(levels(sim$truth), function(label) {
    safe_mean(correct[sim$truth == label])
  }, numeric(1L))
  region <- if (!is.null(sim$region)) sim$region else sim$area
  region_count <- table(droplevels(region))
  sparse_region <- names(which.min(region_count))
  dense_region <- names(which.max(region_count))
  boundary_cut <- stats::quantile(sim$boundary_proximity, 0.20, names = FALSE)
  decision <- attr(pred, "decision")
  data.frame(
    method = method,
    accuracy = mean(correct),
    ari = adjustedRandIndex(pred, sim$truth),
    macro_recall = mean(class_recall, na.rm = TRUE),
    worst_recall = min(class_recall, na.rm = TRUE),
    sparse_region_accuracy = safe_mean(correct[region == sparse_region]),
    dense_region_accuracy = safe_mean(correct[region == dense_region]),
    boundary_accuracy = safe_mean(correct[sim$boundary_proximity <= boundary_cut]),
    interior_accuracy = safe_mean(correct[sim$boundary_proximity > boundary_cut]),
    correction_recall = safe_mean(correct[initially_wrong]),
    damage_rate = safe_mean(!correct[!initially_wrong]),
    changed_precision = safe_mean(correct[changed]),
    changed_fraction = mean(changed),
    unresolved_fraction = if (is.null(decision)) 0 else mean(decision == "unresolved"),
    density_ratio = max(region_count) / max(1, min(region_count)),
    seconds = unname(elapsed["elapsed"]),
    stringsAsFactors = FALSE
  )
}

run_methods <- function(sim, seed) {
  methods <- list(
    "marginSVM" = function() refine_spatial_svm(
      sim$xy, sim$labels, sim$samples,
      control = list(workers = 1L, seed = seed)
    ),
    "marginSVM trust field" = function() refine_spatial_svm(
      sim$xy, sim$labels, sim$samples,
      control = list(experimental_v2 = 1, workers = 1L, seed = seed)
    ),
    "Graph refinement" = function() refine_spatial_clusters(sim$xy, sim$labels, sim$samples),
    "SpaGCN refine" = function() SpatialGraphRefine:::.refine_published_labels(
      sim$xy, sim$labels, sim$samples, method = "spagcn", neighbors = 6L
    ),
    "GraphST refine" = function() SpatialGraphRefine:::.refine_published_labels(
      sim$xy, sim$labels, sim$samples, method = "graphst", neighbors = 50L
    ),
    "C++ kNN vote" = function() refine_spatial_clusters(
      sim$xy, sim$labels, sim$samples,
      control = list(weighted = FALSE, iterations = 1L, consensus = 0.5,
                     preserve = 0, margin = 0, current_support = 1)
    ),
    "Potts-like ICM" = function() refine_spatial_clusters(
      sim$xy, sim$labels, sim$samples,
      control = list(weighted = FALSE, iterations = 8L, consensus = 0.5,
                     preserve = 0.35, margin = 0, current_support = 1)
    )
  )
  output <- list(score_prediction(sim$labels, sim, c(elapsed = 0), "Initial"))
  for (method in names(methods)) {
    timing <- system.time(pred <- methods[[method]]())
    output[[length(output) + 1L]] <- score_prediction(pred, sim, timing, method)
  }
  bind_rows(output)
}

geometry_design <- crossing(
  pattern = patterns,
  density_profile = if (quick) c("uniform", "strong") else density_profiles,
  noise_type = if (quick) c("random", "patch") else noise_types,
  replicate = replicates
) |>
  mutate(
    family = "geometric",
    dimensions = if_else(pattern == "layers3d", 3L, 2L),
    samples = if_else(replicate == max(replicates) & length(replicates) > 1L, 3L, 1L),
    noise = 0.25,
    minority = NA_real_
  )

gradient_design <- crossing(
  dimensions = c(2L, 3L),
  density_profile = if (quick) c("uniform", "strong") else
    c("uniform", "moderate", "strong", "extreme"),
  minority = if (quick) 0.05 else c(0.05, 0.25),
  replicate = replicates
) |>
  mutate(
    family = "gradient",
    pattern = "gradient_A-B-B-C",
    noise_type = "mixture",
    samples = if_else(replicate == max(replicates) & length(replicates) > 1L, 3L, 1L),
    noise = minority
  )

design <- bind_rows(geometry_design, gradient_design) |>
  mutate(scenario = row_number())

run_scenario <- function(index) {
  d <- design[index, ]
  seed <- 1300000L + index
  sim <- if (d$family == "gradient") {
    simulate_gradient_regions(
      n = n_geometry, minority = d$minority, dimensions = d$dimensions,
      samples = d$samples, density_profile = d$density_profile, seed = seed
    )
  } else {
    simulate_spatial_domains(
      n = n_geometry, pattern = d$pattern, dimensions = d$dimensions,
      k = 5L, samples = d$samples, noise = d$noise,
      noise_type = d$noise_type, density_profile = d$density_profile, seed = seed
    )
  }
  metrics <- run_methods(sim, 1400000L + index)
  bind_cols(d[rep(1L, nrow(metrics)), ], metrics)
}

message("Running ", nrow(design), " complete simulation scenarios")
scenario_results <- if (.Platform$OS.type != "windows" && workers > 1L) {
  parallel::mclapply(seq_len(nrow(design)), run_scenario,
                     mc.cores = workers, mc.preschedule = TRUE)
} else {
  lapply(seq_len(nrow(design)), function(i) {
    message("scenario ", i, "/", nrow(design))
    run_scenario(i)
  })
}
metrics <- bind_rows(scenario_results)
metrics$method <- factor(metrics$method, method_levels)
write.csv(metrics, file.path(out_dir, "complete_simulation_metrics.csv"), row.names = FALSE)

summary <- metrics |>
  filter(method != "Initial") |>
  group_by(method) |>
  summarise(
    scenarios = n(), accuracy = mean(accuracy), ari = mean(ari),
    macro_recall = mean(macro_recall), worst_recall = mean(worst_recall),
    sparse_region_accuracy = mean(sparse_region_accuracy),
    boundary_accuracy = mean(boundary_accuracy), correction_recall = mean(correction_recall),
    damage_rate = mean(damage_rate), unresolved_fraction = mean(unresolved_fraction),
    seconds = mean(seconds), .groups = "drop"
  ) |>
  arrange(desc(accuracy))
write.csv(summary, file.path(out_dir, "complete_simulation_summary.csv"), row.names = FALSE)

paired <- metrics |>
  filter(method != "Initial") |>
  select(scenario, method, accuracy) |>
  pivot_wider(names_from = method, values_from = accuracy)
set.seed(1500001L)
comparators <- setdiff(as.character(unique(metrics$method)), c("Initial", "marginSVM"))
comparators <- comparators[!is.na(comparators)]
paired_effects <- bind_rows(lapply(comparators, function(comparator) {
  difference <- paired[["marginSVM"]] - paired[[comparator]]
  bootstrap <- replicate(10000L, mean(sample(difference, replace = TRUE)))
  data.frame(
    comparator = comparator, mean_difference = mean(difference),
    lower_95 = unname(quantile(bootstrap, 0.025)),
    upper_95 = unname(quantile(bootstrap, 0.975)),
    wins = mean(difference > 0), ties = mean(difference == 0)
  )
}))
write.csv(paired_effects, file.path(out_dir, "paired_accuracy_effects.csv"), row.names = FALSE)

density_summary <- metrics |>
  filter(method != "Initial") |>
  group_by(density_profile, noise_type, method) |>
  summarise(accuracy = mean(accuracy), sparse = mean(sparse_region_accuracy),
            damage = mean(damage_rate), density_ratio = mean(density_ratio), .groups = "drop")
write.csv(density_summary, file.path(out_dir, "density_summary.csv"), row.names = FALSE)

p_density <- density_summary |>
  filter(noise_type != "mixture") |>
  mutate(density_profile = factor(density_profile, density_profiles)) |>
  ggplot(aes(density_profile, accuracy, color = method, group = method)) +
  geom_line(linewidth = 0.45) + geom_point(size = 1.2) +
  facet_wrap(~noise_type, nrow = 1L) + coord_cartesian(ylim = c(0.45, 1)) +
  labs(x = "Region concentration profile", y = "Reference accuracy", color = NULL) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
ggsave(file.path(out_dir, "accuracy_by_density_and_error.png"), p_density,
       width = 12.5, height = 4.7, dpi = 240, bg = "white")

p_sparse <- density_summary |>
  filter(noise_type != "mixture") |>
  mutate(density_profile = factor(density_profile, density_profiles)) |>
  ggplot(aes(density_profile, sparse, color = method, group = method)) +
  geom_line(linewidth = 0.45) + geom_point(size = 1.2) +
  facet_wrap(~noise_type, nrow = 1L) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Region concentration profile", y = "Sparsest-region accuracy", color = NULL) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
ggsave(file.path(out_dir, "sparse_region_accuracy.png"), p_sparse,
       width = 12.5, height = 4.7, dpi = 240, bg = "white")

geometry_margin <- metrics |>
  filter(method %in% c("marginSVM", "GraphST refine"), family == "geometric") |>
  select(scenario, pattern, density_profile, noise_type, method, accuracy) |>
  pivot_wider(names_from = method, values_from = accuracy) |>
  mutate(advantage = marginSVM - `GraphST refine`) |>
  group_by(pattern, density_profile) |>
  summarise(advantage = mean(advantage), .groups = "drop") |>
  mutate(density_profile = factor(density_profile, density_profiles))
p_heat <- ggplot(geometry_margin, aes(density_profile, pattern, fill = advantage)) +
  geom_tile() +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C44E52", midpoint = 0) +
  labs(x = "Region concentration profile", y = NULL,
       fill = "marginSVM -\nGraphST") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(out_dir, "geometry_density_advantage.png"), p_heat,
       width = 8.5, height = 5.8, dpi = 240, bg = "white")

gradient_summary <- metrics |>
  filter(family == "gradient", method != "Initial") |>
  group_by(dimensions, density_profile, minority, method) |>
  summarise(accuracy = mean(accuracy), sparse = mean(sparse_region_accuracy), .groups = "drop")
p_gradient <- ggplot(gradient_summary,
                     aes(factor(minority), accuracy, color = method, group = method)) +
  geom_line(linewidth = 0.45) + geom_point(size = 1.2) +
  facet_grid(dimensions ~ density_profile, labeller = label_both) +
  coord_cartesian(ylim = c(0.45, 1)) +
  labs(x = "Within-area minority fraction", y = "Reference accuracy", color = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "gradient_density_accuracy.png"), p_gradient,
       width = 12.5, height = 6.4, dpi = 240, bg = "white")

atlas_patterns <- c(patterns, "gradient_A-B-B-C")
atlas_rows <- list()
atlas_position <- 0L
for (i in seq_along(atlas_patterns)) {
  pattern <- atlas_patterns[i]
  sim <- if (pattern == "gradient_A-B-B-C") {
    simulate_gradient_regions(6000, minority = 0.25, density_profile = "strong", seed = 1600000L + i)
  } else {
    simulate_spatial_domains(
      6000, pattern, dimensions = if (pattern == "layers3d") 3L else 2L,
      noise = 0.25, noise_type = "random", density_profile = "strong",
      seed = 1600000L + i
    )
  }
  pred <- refine_spatial_svm(sim$xy, sim$labels, sim$samples,
                             control = list(workers = 1L, seed = 1700000L + i))
  display <- if (nrow(sim$xy) > 5000L) sample.int(nrow(sim$xy), 5000L) else seq_len(nrow(sim$xy))
  states <- list(Reference = sim$truth, `Corrupted input` = sim$labels, marginSVM = pred)
  for (state in names(states)) {
    atlas_position <- atlas_position + 1L
    atlas_rows[[atlas_position]] <- data.frame(
      x = sim$xy[display, 1L], y = sim$xy[display, 2L],
      pattern = pattern, state = state, label = states[[state]][display]
    )
  }
}
atlas <- bind_rows(atlas_rows)
atlas$state <- factor(atlas$state, c("Reference", "Corrupted input", "marginSVM"))
p_atlas <- ggplot(atlas, aes(x, y, color = label)) +
  geom_point(size = 0.14, alpha = 0.82) + coord_equal() +
  facet_grid(pattern ~ state) + guides(color = "none") +
  labs(x = NULL, y = NULL) + theme_void(base_size = 8) +
  theme(strip.text.y = element_text(angle = 0),
        plot.background = element_rect(fill = "white", color = NA))
ggsave(file.path(out_dir, "complete_geometry_atlas.png"), p_atlas,
       width = 10, height = 22, dpi = 240, bg = "white")

density_examples <- list()
position <- 0L
for (pattern in c("jagged_stripes", "rings", "islands", "thin_layers")) {
  for (profile in density_profiles) {
    sim <- simulate_spatial_domains(
      5000, pattern, noise = 0, density_profile = profile,
      seed = 1800000L + match(pattern, patterns) * 10L + match(profile, density_profiles)
    )
    position <- position + 1L
    density_examples[[position]] <- data.frame(
      x = sim$xy[, 1L], y = sim$xy[, 2L], label = sim$truth,
      pattern = pattern, density_profile = profile
    )
  }
}
density_examples <- bind_rows(density_examples)
density_examples$density_profile <- factor(density_examples$density_profile, density_profiles)
p_examples <- ggplot(density_examples, aes(x, y, color = label)) +
  geom_point(size = 0.17, alpha = 0.8) + coord_equal() +
  facet_grid(pattern ~ density_profile) + guides(color = "none") +
  labs(x = NULL, y = NULL) + theme_void(base_size = 9) +
  theme(strip.text = element_text(face = "bold"),
        plot.background = element_rect(fill = "white", color = NA))
ggsave(file.path(out_dir, "density_profile_examples.png"), p_examples,
       width = 12.5, height = 9.5, dpi = 240, bg = "white")

print(summary)
print(paired_effects)
