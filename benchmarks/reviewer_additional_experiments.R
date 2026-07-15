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

cpp_knn_control <- list(
  weighted = FALSE, iterations = 1L, consensus = 0.5,
  preserve = 0, margin = 0, current_support = 1
)

score_simple <- function(pred, sim, method, seconds) {
  recalls <- vapply(levels(sim$truth), function(level) mean(pred[sim$truth == level] == level), numeric(1))
  data.frame(
    method = method,
    accuracy = mean(pred == sim$truth),
    ari = mclust::adjustedRandIndex(pred, sim$truth),
    smallest_recall = min(recalls),
    damage_rate = mean(pred[sim$labels == sim$truth] != sim$truth[sim$labels == sim$truth]),
    seconds = unname(seconds["elapsed"])
  )
}

run_three <- function(sim) {
  rows <- list(score_simple(sim$labels, sim, "Initial", c(elapsed = 0)))
  timing <- system.time(pred <- refine_spatial_clusters(sim$xy, sim$labels, sim$samples))
  rows[[2L]] <- score_simple(pred, sim, "Spatial graph", timing)
  timing <- system.time(pred <- refine_spatial_clusters(sim$xy, sim$labels, sim$samples, control = cpp_knn_control))
  rows[[3L]] <- score_simple(pred, sim, "C++ kNN vote", timing)
  bind_rows(rows)
}

subset_sim <- function(sim, rows) {
  fields <- c("xy", "labels", "truth", "samples", "boundary_proximity", "corrupted")
  for (field in fields) {
    if (is.matrix(sim[[field]])) sim[[field]] <- sim[[field]][rows, , drop = FALSE]
    else sim[[field]] <- sim[[field]][rows]
  }
  sim
}

message("Running irregular-density and sampling-hole experiment")
density_results <- list()
counter <- 0L
for (pattern in c("jagged_stripes", "islands", "thin_layers")) {
  for (sampling in c("uniform", "gradient", "holes")) {
    for (replicate in seq_len(reps)) {
      counter <- counter + 1L
      sim <- simulate_spatial_domains(30000L, pattern, noise = 0.25, seed = 1000000L + counter)
      x <- sim$xy[, 1L]
      y <- sim$xy[, 2L]
      weights <- switch(
        sampling,
        uniform = rep(1, length(x)),
        gradient = exp(2.2 * x),
        holes = ifelse((x + 0.35)^2 + (y - 0.12)^2 < 0.15^2 |
                       (x - 0.42)^2 + (y + 0.18)^2 < 0.20^2, 0.015, 1)
      )
      rows <- sample.int(nrow(sim$xy), 12000L, prob = weights)
      sim <- subset_sim(sim, rows)
      metrics <- run_three(sim)
      density_results[[counter]] <- bind_cols(
        data.frame(pattern = pattern, sampling = sampling, replicate = replicate), metrics
      )
    }
  }
}
density_metrics <- bind_rows(density_results)
write.csv(density_metrics, file.path(out_dir, "irregular_density_metrics.csv"), row.names = FALSE)

message("Running 3D anisotropy experiment")
anisotropy_results <- list()
counter <- 0L
for (z_scale in c(0.1, 0.25, 1, 4, 10)) {
  for (replicate in seq_len(reps)) {
    counter <- counter + 1L
    sim <- simulate_spatial_domains(15000L, "layers3d", dimensions = 3L, noise = 0.25,
                                    seed = 1100000L + counter)
    physical_xy <- sim$xy
    sim$xy[, 3L] <- sim$xy[, 3L] * z_scale
    timing <- system.time(raw_pred <- refine_spatial_clusters(sim$xy, sim$labels))
    raw <- score_simple(raw_pred, sim, "Uncorrected coordinates", timing)
    sim$xy <- physical_xy
    timing <- system.time(corrected_pred <- refine_spatial_clusters(sim$xy, sim$labels))
    corrected <- score_simple(corrected_pred, sim, "Common physical units", timing)
    anisotropy_results[[counter]] <- bind_cols(
      data.frame(z_scale = z_scale, replicate = replicate), bind_rows(raw, corrected)
    )
  }
}
anisotropy_metrics <- bind_rows(anisotropy_results)
write.csv(anisotropy_metrics, file.path(out_dir, "anisotropy_3d_metrics.csv"), row.names = FALSE)

message("Running duplicated-coordinate experiment")
duplicate_results <- list()
counter <- 0L
for (digits in c(4L, 2L, 1L, 0L)) {
  for (replicate in seq_len(reps)) {
    counter <- counter + 1L
    sim <- simulate_spatial_domains(12000L, "rings", noise = 0.25, seed = 1200000L + counter)
    sim$xy <- round(sim$xy, digits = digits)
    duplicate_fraction <- 1 - nrow(unique(sim$xy)) / nrow(sim$xy)
    metrics <- run_three(sim)
    duplicate_results[[counter]] <- bind_cols(
      data.frame(digits = digits, duplicate_fraction = duplicate_fraction, replicate = replicate), metrics
    )
  }
}
duplicate_metrics <- bind_rows(duplicate_results)
write.csv(duplicate_metrics, file.path(out_dir, "duplicate_coordinate_metrics.csv"), row.names = FALSE)

grid_simulation <- function(side, pattern, noise, seed) {
  set.seed(seed)
  grid <- expand.grid(x = seq(-1, 1, length.out = side), y = seq(-1, 1, length.out = side))
  score <- switch(
    pattern,
    jagged = grid$x + 0.22 * sin(6 * pi * grid$y) + 0.09 * sin(13 * pi * grid$y),
    rings = sqrt(grid$x^2 + grid$y^2) + 0.05 * sin(9 * atan2(grid$y, grid$x)),
    layers = grid$y + 0.28 * sin(2.5 * pi * grid$x)
  )
  truth_id <- as.integer(cut(score, quantile(score, seq(0, 1, length.out = 6)),
                             include.lowest = TRUE, labels = FALSE))
  labels_id <- truth_id
  corrupt <- sample.int(nrow(grid), floor(noise * nrow(grid)))
  replacement <- sample.int(5L, length(corrupt), replace = TRUE)
  replacement[replacement == truth_id[corrupt]] <- replacement[replacement == truth_id[corrupt]] %% 5L + 1L
  labels_id[corrupt] <- replacement
  list(
    xy = as.matrix(grid),
    labels = factor(labels_id, levels = 1:5),
    truth = factor(truth_id, levels = 1:5),
    samples = factor(rep(1L, nrow(grid))),
    side = side
  )
}

modal_filter <- function(labels, side, radius = 1L, iterations = 1L) {
  current <- matrix(as.integer(labels), nrow = side, ncol = side)
  classes <- seq_len(nlevels(labels))
  for (iteration in seq_len(iterations)) {
    votes <- array(0L, dim = c(side, side, length(classes)))
    for (dx in -radius:radius) {
      for (dy in -radius:radius) {
        source_x <- seq.int(max(1L, 1L - dx), min(side, side - dx))
        source_y <- seq.int(max(1L, 1L - dy), min(side, side - dy))
        target_x <- source_x + dx
        target_y <- source_y + dy
        shifted <- current[source_x, source_y, drop = FALSE]
        for (class in classes) votes[target_x, target_y, class] <- votes[target_x, target_y, class] + (shifted == class)
      }
    }
    current <- apply(votes, c(1L, 2L), which.max)
  }
  factor(as.vector(current), levels = classes)
}

message("Running regular-grid morphological comparison")
grid_results <- list()
counter <- 0L
for (pattern in c("jagged", "rings", "layers")) {
  for (noise in c(0.15, 0.25, 0.40)) {
    for (replicate in seq_len(reps)) {
      counter <- counter + 1L
      sim <- grid_simulation(if (quick) 120L else 200L, pattern, noise, 1300000L + counter)
      methods <- list(
        "Initial" = function() sim$labels,
        "Spatial graph" = function() refine_spatial_clusters(sim$xy, sim$labels),
        "SpaGCN refine" = function() SpatialGraphRefine:::.refine_published_labels(
          sim$xy, sim$labels, sim$samples, method = "spagcn", neighbors = 4L
        ),
        "GraphST refine" = function() SpatialGraphRefine:::.refine_published_labels(
          sim$xy, sim$labels, sim$samples, method = "graphst", neighbors = 50L
        ),
        "C++ kNN vote" = function() refine_spatial_clusters(sim$xy, sim$labels, control = cpp_knn_control),
        "3x3 modal filter" = function() modal_filter(sim$labels, sim$side, 1L, 1L),
        "Iterated 3x3 filter" = function() modal_filter(sim$labels, sim$side, 1L, 3L)
      )
      rows <- lapply(names(methods), function(method) {
        timing <- system.time(pred <- methods[[method]]())
        score_simple(pred, sim, method, timing)
      })
      grid_results[[counter]] <- bind_cols(
        data.frame(pattern = pattern, noise = noise, replicate = replicate), bind_rows(rows)
      )
    }
  }
}
grid_metrics <- bind_rows(grid_results)
write.csv(grid_metrics, file.path(out_dir, "regular_grid_metrics.csv"), row.names = FALSE)

message("Computing paired bootstrap intervals")
robust <- read.csv(file.path(out_dir, "robustness_metrics.csv"))
keys <- c("pattern", "noise", "noise_type", "replicate")
wide <- robust |>
  select(all_of(keys), method, accuracy, boundary_010, damage_rate) |>
  pivot_wider(names_from = method, values_from = c(accuracy, boundary_010, damage_rate))

comparators <- c("Initial", "SpaGCN refine", "GraphST refine",
                 "C++ kNN vote", "Potts-like ICM")
set.seed(1400001L)
bootstrap_rows <- list()
counter <- 0L
for (metric in c("accuracy", "boundary_010", "damage_rate")) {
  graph_values <- wide[[paste0(metric, "_Spatial graph")]]
  for (comparator in comparators) {
    difference <- graph_values - wide[[paste0(metric, "_", comparator)]]
    bootstrap <- replicate(5000L, mean(sample(difference, replace = TRUE), na.rm = TRUE))
    counter <- counter + 1L
    bootstrap_rows[[counter]] <- data.frame(
      metric = metric, comparator = comparator,
      mean_difference = mean(difference, na.rm = TRUE),
      lower_95 = quantile(bootstrap, 0.025, na.rm = TRUE),
      upper_95 = quantile(bootstrap, 0.975, na.rm = TRUE)
    )
  }
}
bootstrap_metrics <- bind_rows(bootstrap_rows)
write.csv(bootstrap_metrics, file.path(out_dir, "paired_bootstrap_intervals.csv"), row.names = FALSE)

p_density <- density_metrics |>
  group_by(pattern, sampling, method) |>
  summarise(accuracy = mean(accuracy), .groups = "drop") |>
  ggplot(aes(sampling, accuracy, color = method, group = method)) +
  geom_point() + geom_line() + facet_wrap(~pattern) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Sampling design", y = "Accuracy", color = "Method") + theme_bw(base_size = 10)
ggsave(file.path(out_dir, "irregular_density_accuracy.png"), p_density, width = 9, height = 4, dpi = 220)

p_anisotropy <- anisotropy_metrics |>
  group_by(z_scale, method) |>
  summarise(accuracy = mean(accuracy), .groups = "drop") |>
  ggplot(aes(z_scale, accuracy, color = method)) +
  geom_point() + geom_line() + scale_x_log10() + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Numerical z-axis scale", y = "3D accuracy", color = NULL) + theme_bw(base_size = 10)
ggsave(file.path(out_dir, "anisotropy_3d_accuracy.png"), p_anisotropy, width = 7, height = 4, dpi = 220)

p_grid <- grid_metrics |>
  group_by(pattern, noise, method) |>
  summarise(accuracy = mean(accuracy), .groups = "drop") |>
  ggplot(aes(noise, accuracy, color = method)) +
  geom_point() + geom_line() + facet_wrap(~pattern) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Corrupted-label fraction", y = "Regular-grid accuracy", color = "Method") +
  theme_bw(base_size = 10)
ggsave(file.path(out_dir, "regular_grid_comparison.png"), p_grid, width = 9, height = 4, dpi = 220)

message("Additional reviewer-response experiments complete")
