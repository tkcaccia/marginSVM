#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SpatialGraphRefine)
  library(FNN)
  library(mclust)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

object_path <- Sys.getenv(
  "SPATIAL_REFINE_CRC_RDS",
  "/Users/stefano/Documents/wsitools/Data/VisiumHD/Colorectal/SpatialPolygons_2000genes_counts_only_with_tissue_annotations.rds"
)
out_dir <- file.path("benchmarks", "results", "visiumhd_colorectal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workers <- max(1L, as.integer(Sys.getenv("SPATIAL_REFINE_BENCHMARK_WORKERS", "4")))
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank()))

message("Loading annotated VisiumHD colorectal section")
object <- readRDS(object_path)
xy_all <- as.matrix(GetTissueCoordinates(object)[, c("x", "y")])
labels_all <- object@meta.data$wsi_annotation_name
if (!identical(rownames(xy_all), rownames(object@meta.data))) {
  stop("Coordinate and metadata row orders do not match.")
}
keep <- !is.na(labels_all) & is.finite(xy_all[, 1L]) & is.finite(xy_all[, 2L])
xy <- xy_all[keep, , drop = FALSE]
truth <- droplevels(factor(labels_all[keep]))
truth_code <- as.integer(truth)
n <- nrow(xy)
n_classes <- nlevels(truth)

write.csv(
  data.frame(
    observations_total = nrow(xy_all), observations_annotated = n,
    observations_excluded = sum(!keep), classes = n_classes,
    duplicated_coordinates = sum(duplicated(xy))
  ),
  file.path(out_dir, "dataset_summary.csv"), row.names = FALSE
)
write.csv(
  data.frame(annotation = names(table(truth)), observations = as.integer(table(truth))),
  file.path(out_dir, "annotation_counts.csv"), row.names = FALSE
)

message("Building annotation-boundary reference graph")
nn8 <- FNN::get.knn(xy, k = 8L)$nn.index
neighbor_codes <- matrix(truth_code[nn8], nrow = nrow(nn8))
boundary <- rowSums(neighbor_codes != truth_code) > 0L
adjacent_code <- truth_code
for (column in seq_len(ncol(nn8))) {
  candidate <- neighbor_codes[, column]
  replace <- adjacent_code == truth_code & candidate != truth_code
  adjacent_code[replace] <- candidate[replace]
}
adjacent_code[adjacent_code == truth_code] <- (truth_code[adjacent_code == truth_code] %% n_classes) + 1L

unit_xy <- sweep(xy, 2L, apply(xy, 2L, min), "-")
unit_xy <- sweep(unit_xy, 2L, pmax(apply(unit_xy, 2L, max), .Machine$double.eps), "/")

different_random_codes <- function(indices) {
  replacement <- sample.int(n_classes - 1L, length(indices), replace = TRUE)
  replacement + (replacement >= truth_code[indices])
}

corrupt_labels <- function(fraction, mechanism, seed) {
  if (fraction <= 0) return(truth)
  set.seed(seed)
  n_noise <- min(n - 1L, as.integer(round(fraction * n)))
  if (mechanism == "random") {
    indices <- sample.int(n, n_noise)
    replacement <- different_random_codes(indices)
  } else if (mechanism == "boundary") {
    boundary_rows <- which(boundary)
    if (length(boundary_rows) >= n_noise) {
      indices <- sample(boundary_rows, n_noise)
    } else {
      indices <- c(boundary_rows, sample(setdiff(seq_len(n), boundary_rows), n_noise - length(boundary_rows)))
    }
    replacement <- adjacent_code[indices]
  } else if (mechanism == "patch") {
    centers <- sample.int(n, 12L)
    distance_matrix <- vapply(
      centers,
      function(center) rowSums((unit_xy - unit_xy[center, ])^2),
      numeric(n)
    )
    owner <- max.col(-distance_matrix, ties.method = "first")
    distance <- distance_matrix[cbind(seq_len(n), owner)]
    indices <- order(distance)[seq_len(n_noise)]
    center_targets <- adjacent_code[centers]
    replacement <- center_targets[owner[indices]]
    same <- replacement == truth_code[indices]
    replacement[same] <- adjacent_code[indices[same]]
  } else if (mechanism == "region") {
    center <- sample.int(n, 1L)
    distance <- rowSums((unit_xy - unit_xy[center, ])^2)
    indices <- order(distance)[seq_len(n_noise)]
    replacement <- rep.int(adjacent_code[center], n_noise)
    same <- replacement == truth_code[indices]
    replacement[same] <- adjacent_code[indices[same]]
  } else {
    stop("Unknown corruption mechanism: ", mechanism)
  }
  output <- truth
  output[indices] <- levels(truth)[replacement]
  attr(output, "corrupted") <- indices
  output
}

safe_mean <- function(x) if (length(x)) mean(x) else NA_real_

score_prediction <- function(pred, initial, method, elapsed) {
  correct <- pred == truth
  initially_wrong <- initial != truth
  initially_correct <- !initially_wrong
  changed <- pred != initial
  recalls <- vapply(
    levels(truth),
    function(annotation) safe_mean(pred[truth == annotation] == annotation),
    numeric(1L)
  )
  data.frame(
    method = method,
    accuracy = mean(correct),
    ari = mclust::adjustedRandIndex(pred, truth),
    macro_recall = mean(recalls),
    worst_recall = min(recalls),
    boundary_accuracy = safe_mean(correct[boundary]),
    interior_accuracy = safe_mean(correct[!boundary]),
    changed_precision = safe_mean(correct[changed]),
    correction_recall = safe_mean(correct[initially_wrong]),
    damage_rate = safe_mean(!correct[initially_correct]),
    changed_fraction = mean(changed),
    seconds = unname(elapsed),
    stringsAsFactors = FALSE
  )
}

control_from_row <- function(row) {
  if (is.na(row$neighbors)) return(list())
  list(neighbors = as.integer(row$neighbors), current_support = row$current_support)
}

candidate_grid <- bind_rows(
  data.frame(candidate = "Legacy default", neighbors = 28L, current_support = 0.25),
  crossing(neighbors = c(11L, 15L, 21L, 28L), current_support = c(0.10, 0.15, 0.20)) |>
    mutate(candidate = sprintf("k%d_support%.2f", neighbors, current_support)) |>
    select(candidate, neighbors, current_support)
)

message("Evaluating clean-label preservation for candidate profiles")
clean_candidates <- bind_rows(lapply(seq_len(nrow(candidate_grid)), function(i) {
  row <- candidate_grid[i, ]
  timing <- system.time(pred <- refine_spatial_clusters(xy, truth, control = control_from_row(row)))
  data.frame(
    row, clean_accuracy = mean(pred == truth), clean_damage = mean(pred != truth),
    seconds = unname(timing["elapsed"])
  )
}))

development_design <- crossing(
  noise = c(0.10, 0.25, 0.40),
  mechanism = c("random", "boundary", "patch"),
  replicate = 1:3
)

message("Running development-only parameter selection")
run_development <- function(index) {
  scenario <- development_design[index, ]
  initial <- corrupt_labels(
    scenario$noise, scenario$mechanism,
    seed = 810000L + index
  )
  rows <- lapply(seq_len(nrow(candidate_grid)), function(i) {
    candidate <- candidate_grid[i, ]
    timing <- system.time(pred <- refine_spatial_clusters(
      xy, initial, control = control_from_row(candidate)
    ))
    data.frame(
      scenario, candidate = candidate$candidate,
      accuracy = mean(pred == truth),
      damage_rate = mean(pred[initial == truth] != truth[initial == truth]),
      correction_recall = mean(pred[initial != truth] == truth[initial != truth]),
      seconds = unname(timing["elapsed"])
    )
  })
  bind_rows(rows)
}

development_results <- if (.Platform$OS.type != "windows" && workers > 1L) {
  parallel::mclapply(seq_len(nrow(development_design)), run_development,
                     mc.cores = workers, mc.preschedule = TRUE)
} else {
  lapply(seq_len(nrow(development_design)), run_development)
}
development_metrics <- bind_rows(development_results)
development_summary <- development_metrics |>
  group_by(candidate) |>
  summarise(
    development_accuracy = mean(accuracy),
    development_damage = mean(damage_rate),
    development_correction = mean(correction_recall),
    .groups = "drop"
  ) |>
  left_join(clean_candidates, by = "candidate") |>
  arrange(desc(development_accuracy), clean_damage)
eligible <- development_summary |> filter(clean_damage <= 0.005)
if (!nrow(eligible)) stop("No candidate met the pre-specified 0.5% clean-damage constraint.")
selected_name <- eligible$candidate[1L]
selected_row <- candidate_grid[candidate_grid$candidate == selected_name, , drop = FALSE]
selected_control <- control_from_row(selected_row)
write.csv(development_metrics, file.path(out_dir, "development_metrics.csv"), row.names = FALSE)
write.csv(development_summary, file.path(out_dir, "development_selection.csv"), row.names = FALSE)
write.csv(selected_row, file.path(out_dir, "selected_profile.csv"), row.names = FALSE)
message("Selected profile: ", selected_name)

methods <- list(
  "Legacy default" = function(initial) refine_spatial_clusters(
    xy, initial, control = list(neighbors = 28L, current_support = 0.25)
  ),
  "Adaptive default" = function(initial) refine_spatial_clusters(xy, initial),
  "Selected profile" = function(initial) refine_spatial_clusters(xy, initial, control = selected_control),
  "SpaGCN refine" = function(initial) SpatialGraphRefine:::.refine_published_labels(
    xy, initial, factor(rep.int(1L, n)), method = "spagcn", neighbors = 6L
  ),
  "GraphST refine" = function(initial) SpatialGraphRefine:::.refine_published_labels(
    xy, initial, factor(rep.int(1L, n)), method = "graphst", neighbors = 50L
  ),
  "C++ kNN vote" = function(initial) refine_spatial_clusters(
    xy, initial,
    control = list(weighted = FALSE, iterations = 1L, consensus = 0.5,
                   preserve = 0, margin = 0, current_support = 1)
  ),
  "Potts-like ICM" = function(initial) refine_spatial_clusters(
    xy, initial,
    control = list(weighted = FALSE, iterations = 8L, consensus = 0.5,
                   preserve = 0.35, margin = 0, current_support = 1)
  )
)

test_design <- crossing(
  noise = c(0.05, 0.10, 0.15, 0.25, 0.40),
  mechanism = c("random", "boundary", "patch", "region"),
  replicate = 4:8
)

message("Running held-out real-tissue corruption benchmark")
run_test <- function(index) {
  scenario <- test_design[index, ]
  initial <- corrupt_labels(
    scenario$noise, scenario$mechanism,
    seed = 910000L + index
  )
  initial_row <- score_prediction(initial, initial, "Initial", 0)
  rows <- lapply(names(methods), function(method) {
    timing <- system.time(pred <- methods[[method]](initial))
    score_prediction(pred, initial, method, timing["elapsed"])
  })
  bind_cols(scenario[rep(1L, length(rows) + 1L), ], bind_rows(c(list(initial_row), rows)))
}

test_results <- if (.Platform$OS.type != "windows" && workers > 1L) {
  parallel::mclapply(seq_len(nrow(test_design)), run_test,
                     mc.cores = workers, mc.preschedule = TRUE)
} else {
  lapply(seq_len(nrow(test_design)), run_test)
}
test_metrics <- bind_rows(test_results)

message("Running clean-label negative control")
clean_rows <- bind_rows(c(
  list(score_prediction(truth, truth, "Initial", 0)),
  lapply(names(methods), function(method) {
    timing <- system.time(pred <- methods[[method]](truth))
    score_prediction(pred, truth, method, timing["elapsed"])
  })
)) |>
  mutate(noise = 0, mechanism = "clean", replicate = 0L, .before = 1L)
test_metrics <- bind_rows(test_metrics, clean_rows)
write.csv(test_metrics, file.path(out_dir, "heldout_metrics.csv"), row.names = FALSE)

summary_metrics <- test_metrics |>
  filter(mechanism != "clean") |>
  group_by(method) |>
  summarise(
    accuracy = mean(accuracy), ari = mean(ari), macro_recall = mean(macro_recall),
    worst_recall = mean(worst_recall), boundary_accuracy = mean(boundary_accuracy),
    damage_rate = mean(damage_rate), correction_recall = mean(correction_recall),
    seconds = mean(seconds), .groups = "drop"
  ) |>
  arrange(desc(accuracy))
write.csv(summary_metrics, file.path(out_dir, "heldout_summary.csv"), row.names = FALSE)

paired <- test_metrics |>
  filter(method %in% c("Legacy default", "Adaptive default"), mechanism != "clean") |>
  select(noise, mechanism, replicate, method, accuracy, boundary_accuracy, damage_rate) |>
  pivot_wider(names_from = method, values_from = c(accuracy, boundary_accuracy, damage_rate))
set.seed(920001L)
paired_intervals <- bind_rows(lapply(c("accuracy", "boundary_accuracy", "damage_rate"), function(metric) {
  difference <- paired[[paste0(metric, "_Adaptive default")]] -
    paired[[paste0(metric, "_Legacy default")]]
  bootstrap <- replicate(10000L, mean(sample(difference, replace = TRUE)))
  data.frame(
    metric = metric, mean_difference = mean(difference),
    lower_95 = unname(quantile(bootstrap, 0.025)),
    upper_95 = unname(quantile(bootstrap, 0.975))
  )
}))
write.csv(paired_intervals, file.path(out_dir, "adaptive_vs_default_paired.csv"), row.names = FALSE)

plot_summary <- test_metrics |>
  filter(method %in% c("Initial", "Legacy default", "Adaptive default",
                       "SpaGCN refine", "GraphST refine"), mechanism != "clean") |>
  group_by(noise, mechanism, method) |>
  summarise(accuracy = mean(accuracy), damage_rate = mean(damage_rate), .groups = "drop")

p_accuracy <- ggplot(plot_summary, aes(noise, accuracy, color = method)) +
  geom_line() + geom_point(size = 1.5) + facet_wrap(~mechanism) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(x = "Mislabeled fraction", y = "Held-out annotation accuracy", color = "Method") +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "heldout_accuracy.png"), p_accuracy,
       width = 10, height = 5.7, dpi = 220)

p_damage <- ggplot(
  plot_summary |> filter(method != "Initial"),
  aes(noise, damage_rate, color = method)
) +
  geom_line() + geom_point(size = 1.5) + facet_wrap(~mechanism, scales = "free_y") +
  labs(x = "Mislabeled fraction", y = "Damage among initially correct labels", color = "Method") +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "heldout_damage.png"), p_damage,
       width = 10, height = 5.7, dpi = 220)

message("Generating spatial example and class-level recall")
example_initial <- corrupt_labels(0.25, "patch", 930001L)
example_legacy <- refine_spatial_clusters(
  xy, example_initial, control = list(neighbors = 28L, current_support = 0.25)
)
example_adaptive <- refine_spatial_clusters(xy, example_initial)
set.seed(930002L)
display_rows <- sort(sample.int(n, min(n, 70000L)))
map_data <- bind_rows(
  data.frame(xy[display_rows, ], panel = "Reference", label = truth[display_rows]),
  data.frame(xy[display_rows, ], panel = "25% patch errors", label = example_initial[display_rows]),
  data.frame(xy[display_rows, ], panel = "Legacy default", label = example_legacy[display_rows]),
  data.frame(xy[display_rows, ], panel = "Adaptive default", label = example_adaptive[display_rows])
)
map_data$panel <- factor(map_data$panel, c("Reference", "25% patch errors", "Legacy default", "Adaptive default"))
p_map <- ggplot(map_data, aes(x, y, color = label)) +
  geom_point(size = 0.18, alpha = 0.8) + facet_wrap(~panel, ncol = 2) +
  coord_equal() + scale_y_reverse() +
  labs(x = NULL, y = NULL, color = "Annotation") +
  theme_void(base_size = 9) +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )
ggsave(file.path(out_dir, "spatial_patch_example.png"), p_map,
       width = 11, height = 8, dpi = 240, bg = "white")

class_recall <- bind_rows(lapply(
  list("Initial" = example_initial, "Legacy default" = example_legacy,
       "Adaptive default" = example_adaptive),
  function(pred) data.frame(
    annotation = levels(truth),
    recall = vapply(levels(truth), function(level) mean(pred[truth == level] == level), numeric(1L))
  )
), .id = "method")
write.csv(class_recall, file.path(out_dir, "patch_example_class_recall.csv"), row.names = FALSE)
p_recall <- ggplot(class_recall, aes(reorder(annotation, recall), recall, color = method)) +
  geom_point(position = position_dodge(width = 0.55), size = 1.8) +
  coord_flip(ylim = c(0, 1)) +
  labs(x = NULL, y = "Recall in 25% patch example", color = "Method") +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "patch_example_class_recall.png"), p_recall,
       width = 8, height = 6.5, dpi = 220)

message("VisiumHD colorectal benchmark complete: ", out_dir)
