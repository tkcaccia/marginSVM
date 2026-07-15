#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialGraphRefine)
  library(e1071)
  library(mclust)
  library(ggplot2)
  library(dplyr)
})

out_dir <- file.path("benchmarks", "results", "reviewer_response")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
quick <- identical(Sys.getenv("SPATIAL_REFINE_QUICK"), "1")

tile_exact_svm <- function(xy, labels, tiles, gamma, cost) {
  breaks <- lapply(seq_len(ncol(xy)), function(j) seq(min(xy[, j]), max(xy[, j]), length.out = tiles + 1L))
  tile_id <- rep(0L, nrow(xy))
  multiplier <- 1L
  for (j in seq_len(ncol(xy))) {
    position <- pmin(tiles, findInterval(xy[, j], breaks[[j]], all.inside = TRUE))
    tile_id <- tile_id + (position - 1L) * multiplier
    multiplier <- multiplier * tiles
  }
  out <- labels
  for (tile in unique(tile_id)) {
    rows <- which(tile_id == tile)
    if (length(unique(labels[rows])) < 2L) next
    model <- e1071::svm(
      xy[rows, , drop = FALSE], labels[rows], kernel = "radial",
      gamma = gamma, cost = cost, scale = FALSE
    )
    out[rows] <- predict(model, xy[rows, , drop = FALSE])
  }
  factor(out, levels = levels(labels))
}

validation_patterns <- c("jagged_stripes", "rings", "islands")
gamma_grid <- c(0.1, 0.5, 1, 2)
cost_grid <- c(0.5, 1, 5)
tile_grid <- c(4L, 6L, 8L)
n_validation <- if (quick) 2500L else 4000L
validation_reps <- if (quick) 1L else 3L

global_tuning <- expand.grid(gamma = gamma_grid, cost = cost_grid)
global_tuning$accuracy <- 0
tile_tuning <- expand.grid(gamma = gamma_grid, cost = cost_grid, tiles = tile_grid)
tile_tuning$accuracy <- 0
rff_tuning <- expand.grid(gamma = gamma_grid, tiles = tile_grid)
rff_tuning$accuracy <- 0

validation_data <- list()
counter <- 0L
for (pattern in validation_patterns) {
  for (replicate in seq_len(validation_reps)) {
    counter <- counter + 1L
    validation_data[[counter]] <- simulate_spatial_domains(
      n_validation, pattern, noise = 0.25, noise_type = "random", seed = 900000L + counter
    )
  }
}

message("Tuning SVM comparators on validation simulations")
for (i in seq_len(nrow(global_tuning))) {
  scores <- vapply(validation_data, function(sim) {
    model <- e1071::svm(sim$xy, sim$labels, kernel = "radial",
                        gamma = global_tuning$gamma[i], cost = global_tuning$cost[i], scale = FALSE)
    mean(predict(model, sim$xy) == sim$truth)
  }, numeric(1))
  global_tuning$accuracy[i] <- mean(scores)
}
for (i in seq_len(nrow(tile_tuning))) {
  scores <- vapply(validation_data, function(sim) {
    pred <- tile_exact_svm(sim$xy, sim$labels, tile_tuning$tiles[i],
                           tile_tuning$gamma[i], tile_tuning$cost[i])
    mean(pred == sim$truth)
  }, numeric(1))
  tile_tuning$accuracy[i] <- mean(scores)
}
for (i in seq_len(nrow(rff_tuning))) {
  scores <- vapply(validation_data, function(sim) {
    pred <- refine_spatial_svm(
      sim$xy, sim$labels, tiles = rep(rff_tuning$tiles[i], 2L), gamma = rff_tuning$gamma[i],
      n_features = 96L, epochs = 6L, seed = 42L
    )
    mean(pred == sim$truth)
  }, numeric(1))
  rff_tuning$accuracy[i] <- mean(scores)
}

best_global <- global_tuning[which.max(global_tuning$accuracy), ]
best_tile <- tile_tuning[which.max(tile_tuning$accuracy), ]
best_rff <- rff_tuning[which.max(rff_tuning$accuracy), ]
write.csv(global_tuning, file.path(out_dir, "global_svm_validation.csv"), row.names = FALSE)
write.csv(tile_tuning, file.path(out_dir, "tile_svm_validation.csv"), row.names = FALSE)
write.csv(rff_tuning, file.path(out_dir, "rff_svm_validation.csv"), row.names = FALSE)

message("Testing frozen SVM settings on held-out simulations")
test_patterns <- c("wavy_layers", "spiral", "branching", "lobes", "thin_layers", "disconnected")
test_reps <- if (quick) 2L else 10L
results <- list()
counter <- 0L
for (pattern in test_patterns) {
  for (replicate in seq_len(test_reps)) {
    sim <- simulate_spatial_domains(n_validation, pattern, noise = 0.25, noise_type = "random",
                                    seed = 910000L + counter)
    methods <- list(
      "Spatial graph" = function() refine_spatial_clusters(sim$xy, sim$labels),
      "C++ kNN vote" = function() refine_spatial_clusters(
        sim$xy, sim$labels,
        control = list(weighted = FALSE, iterations = 1L, consensus = 0.5,
                       preserve = 0, margin = 0, current_support = 1)
      ),
      "Global exact SVM" = function() {
        model <- e1071::svm(sim$xy, sim$labels, kernel = "radial",
                            gamma = best_global$gamma, cost = best_global$cost, scale = FALSE)
        factor(predict(model, sim$xy), levels = levels(sim$labels))
      },
      "Tile exact SVM" = function() tile_exact_svm(
        sim$xy, sim$labels, best_tile$tiles, best_tile$gamma, best_tile$cost
      ),
      "Tile RFF-SVM" = function() refine_spatial_svm(
        sim$xy, sim$labels, tiles = rep(best_rff$tiles, 2L), gamma = best_rff$gamma,
        n_features = 96L, epochs = 6L, seed = 42L
      )
    )
    for (method in names(methods)) {
      counter <- counter + 1L
      timing <- system.time(pred <- methods[[method]]())
      results[[counter]] <- data.frame(
        pattern = pattern, replicate = replicate, method = method,
        accuracy = mean(pred == sim$truth),
        ari = mclust::adjustedRandIndex(pred, sim$truth),
        seconds = unname(timing["elapsed"])
      )
    }
  }
}
metrics <- bind_rows(results)
write.csv(metrics, file.path(out_dir, "svm_comparator_test.csv"), row.names = FALSE)
write.csv(bind_rows(
  transform(best_global, method = "Global exact SVM"),
  transform(best_tile, method = "Tile exact SVM"),
  transform(best_rff, method = "Tile RFF-SVM")
), file.path(out_dir, "svm_selected_parameters.csv"), row.names = FALSE)

summary <- metrics |>
  group_by(method) |>
  summarise(accuracy = mean(accuracy), ari = mean(ari), seconds = mean(seconds), .groups = "drop")
write.csv(summary, file.path(out_dir, "svm_comparator_summary.csv"), row.names = FALSE)

p <- ggplot(metrics, aes(method, accuracy, color = method)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, size = 0.8) +
  coord_cartesian(ylim = c(0, 1)) + facet_wrap(~pattern) +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none") +
  labs(x = NULL, y = "Held-out assignment accuracy")
ggsave(file.path(out_dir, "svm_comparator_accuracy.png"), p, width = 11, height = 6, dpi = 220)
