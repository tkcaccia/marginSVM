#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialGraphRefine)
  library(ggplot2)
})

output_dir <- file.path("benchmarks", "results", "tiled_execution")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

time_refinement <- function(sim, execution, repeat_id) {
  gc()
  elapsed <- system.time({
    prediction <- refine_spatial_clusters(
      sim$xy, sim$labels, sim$samples, execution = execution
    )
  })[["elapsed"]]
  list(prediction = prediction, elapsed = elapsed, repeat_id = repeat_id)
}

set.seed(20260714)
scenarios <- list(
  `2D` = simulate_spatial_domains(
    150000, "jagged_stripes", noise = 0.25, samples = 2, seed = 801
  ),
  `3D` = simulate_spatial_domains(
    150000, "layers3d", dimensions = 3, noise = 0.25, samples = 2, seed = 802
  )
)

runtime_rows <- list()
quality_rows <- list()
row_id <- 0L
for (dimension in names(scenarios)) {
  sim <- scenarios[[dimension]]
  full <- refine_spatial_clusters(sim$xy, sim$labels, sim$samples)
  tile_shape <- if (dimension == "2D") c(3, 3) else c(3, 2, 2)
  configurations <- list(
    `Untiled CPU` = NULL,
    `Tiled, 1 CPU` = list(tiles = tile_shape, overlap = 0.15, workers = 1),
    `Tiled, 2 CPUs` = list(tiles = tile_shape, overlap = 0.15, workers = 2),
    `Tiled, 4 CPUs` = list(tiles = tile_shape, overlap = 0.15, workers = 4)
  )
  for (name in names(configurations)) {
    for (repeat_id in seq_len(3L)) {
      measured <- time_refinement(sim, configurations[[name]], repeat_id)
      row_id <- row_id + 1L
      runtime_rows[[row_id]] <- data.frame(
        dimensions = dimension, method = name, replicate = repeat_id,
        seconds = measured$elapsed,
        accuracy = mean(measured$prediction == sim$truth),
        agreement_full = mean(measured$prediction == full),
        workers = attr(measured$prediction, "workers"),
        tile_count = if (is.null(attr(measured$prediction, "tiles"))) 1L else nrow(attr(measured$prediction, "tiles"))
      )
    }
  }

  for (overlap in c(0, 0.05, 0.15, 0.30)) {
    prediction <- refine_spatial_clusters(
      sim$xy, sim$labels, sim$samples,
      execution = list(tiles = tile_shape, overlap = overlap, workers = 1)
    )
    quality_rows[[length(quality_rows) + 1L]] <- data.frame(
      dimensions = dimension, overlap = overlap,
      accuracy = mean(prediction == sim$truth),
      agreement_full = mean(prediction == full),
      halo_multiplier = sum(attr(prediction, "tiles")$halo_n) / nrow(sim$xy)
    )
  }
}

runtime <- do.call(rbind, runtime_rows)
quality <- do.call(rbind, quality_rows)
write.csv(runtime, file.path(output_dir, "tiled_runtime_metrics.csv"), row.names = FALSE)
write.csv(quality, file.path(output_dir, "overlap_quality_metrics.csv"), row.names = FALSE)

runtime_summary <- aggregate(
  cbind(seconds, accuracy, agreement_full) ~ dimensions + method + workers + tile_count,
  runtime, mean
)
write.csv(runtime_summary, file.path(output_dir, "tiled_runtime_summary.csv"), row.names = FALSE)

method_order <- c("Untiled CPU", "Tiled, 1 CPU", "Tiled, 2 CPUs", "Tiled, 4 CPUs")
runtime$method <- factor(runtime$method, levels = method_order)
p_runtime <- ggplot(runtime, aes(method, seconds, color = dimensions)) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.55, linewidth = 0.45) +
  labs(x = NULL, y = "Elapsed seconds", color = "Coordinates") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top")
ggsave(file.path(output_dir, "tiled_multicore_runtime.png"), p_runtime, width = 7.2, height = 4.6, dpi = 180)

p_overlap <- ggplot(quality, aes(overlap, agreement_full, color = dimensions)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  scale_y_continuous(
    limits = c(0.98, 1),
    labels = function(x) paste0(formatC(100 * x, format = "f", digits = 1), "%")
  ) +
  labs(x = "Halo overlap", y = "Agreement with untiled refinement", color = "Coordinates") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
ggsave(file.path(output_dir, "tile_overlap_agreement.png"), p_overlap, width = 7.2, height = 4.6, dpi = 180)

print(runtime_summary)
print(quality)
