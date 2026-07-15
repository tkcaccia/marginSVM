#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialGraphRefine)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(FNN)
  library(mclust)
})

out_dir <- file.path("benchmarks", "results", "publication")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank()))

quick <- identical(Sys.getenv("SPATIAL_REFINE_QUICK"), "1")
patterns <- c(
  "jagged_stripes", "wavy_layers", "rings", "spiral",
  "branching", "lobes", "islands", "disconnected"
)
design <- tidyr::crossing(
  pattern = patterns,
  noise = if (quick) 0.25 else c(0.05, 0.15, 0.25, 0.40),
  noise_type = if (quick) "random" else c("random", "boundary", "patch"),
  replicate = seq_len(if (quick) 3L else 20L)
) |>
  mutate(n = if (quick) 12000L else 20000L, k = 5L)

majority_refine <- function(xy, labels, k = 20L) {
  nn <- FNN::get.knn(xy, k = min(k, nrow(xy) - 1L))$nn.index
  factor(
    apply(nn, 1L, function(i) names(which.max(table(labels[i])))),
    levels = levels(labels)
  )
}

boundary_indicator <- function(xy, truth, k = 8L) {
  nn <- FNN::get.knn(xy, k = min(k, nrow(xy) - 1L))$nn.index
  vapply(seq_len(nrow(nn)), function(i) any(truth[nn[i, ]] != truth[i]), logical(1))
}

score_prediction <- function(pred, truth, boundary) {
  recalls <- vapply(levels(truth), function(level) mean(pred[truth == level] == level), numeric(1))
  c(
    accuracy = mean(pred == truth),
    ari = mclust::adjustedRandIndex(pred, truth),
    macro_recall = mean(recalls),
    rare_recall = min(recalls),
    boundary_accuracy = mean(pred[boundary] == truth[boundary]),
    interior_accuracy = mean(pred[!boundary] == truth[!boundary])
  )
}

run_one <- function(method, expression, truth, boundary) {
  gc(FALSE)
  timing <- system.time(pred <- force(expression))
  as.data.frame(as.list(c(method = method, seconds = unname(timing["elapsed"]),
                          score_prediction(pred, truth, boundary))))
}

results <- vector("list", nrow(design))
for (i in seq_len(nrow(design))) {
  scenario <- design[i, ]
  message("[", i, "/", nrow(design), "] ", scenario$pattern,
          ", noise=", scenario$noise, ", ", scenario$noise_type,
          ", replicate=", scenario$replicate)
  sim <- simulate_spatial_domains(
    n = scenario$n,
    pattern = scenario$pattern,
    k = scenario$k,
    noise = scenario$noise,
    noise_type = scenario$noise_type,
    seed = 10000L + i
  )
  boundary <- boundary_indicator(sim$xy, sim$truth)
  rows <- bind_rows(
    run_one("Initial", sim$labels, sim$truth, boundary),
    run_one("Spatial graph", refine_spatial_clusters(sim$xy, sim$labels), sim$truth, boundary),
    run_one("kNN majority", majority_refine(sim$xy, sim$labels), sim$truth, boundary),
    run_one(
      "Tile RFF-SVM",
      refine_spatial_svm(
        sim$xy, sim$labels, tiles = c(8L, 8L), gamma = 0.65,
        n_features = 96L, epochs = 6L, seed = 42L
      ),
      sim$truth,
      boundary
    )
  )
  results[[i]] <- bind_cols(scenario[rep(1L, nrow(rows)), ], rows)
}

metrics <- bind_rows(results) |>
  mutate(across(c(seconds, accuracy, ari, macro_recall, rare_recall,
                  boundary_accuracy, interior_accuracy), as.numeric))
write.csv(metrics, file.path(out_dir, "geometric_benchmark_metrics.csv"), row.names = FALSE)

summary_metrics <- metrics |>
  group_by(pattern, noise, noise_type, method) |>
  summarise(
    across(c(seconds, accuracy, ari, macro_recall, rare_recall,
             boundary_accuracy, interior_accuracy),
           list(mean = mean, sd = sd), .names = "{.col}_{.fn}"),
    .groups = "drop"
  )
write.csv(summary_metrics, file.path(out_dir, "geometric_benchmark_summary.csv"), row.names = FALSE)

p_accuracy <- ggplot(summary_metrics, aes(pattern, accuracy_mean, color = method, group = method)) +
  geom_point(position = position_dodge(width = 0.35), size = 2) +
  geom_errorbar(
    aes(ymin = accuracy_mean - accuracy_sd, ymax = accuracy_mean + accuracy_sd),
    position = position_dodge(width = 0.35), width = 0.15
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(noise_type ~ noise, labeller = label_both) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = NULL, y = "Assignment accuracy", color = "Method")
ggsave(file.path(out_dir, "geometric_accuracy.png"), p_accuracy, width = 12, height = 6, dpi = 220)

p_boundary <- summary_metrics |>
  select(pattern, noise, noise_type, method, boundary_accuracy_mean, interior_accuracy_mean) |>
  pivot_longer(ends_with("accuracy_mean"), names_to = "region", values_to = "accuracy") |>
  mutate(region = recode(region, boundary_accuracy_mean = "Boundary", interior_accuracy_mean = "Interior")) |>
  ggplot(aes(pattern, accuracy, color = method, group = method)) +
  geom_point(position = position_dodge(width = 0.35), size = 1.8) +
  facet_grid(region ~ noise_type) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom") +
  labs(x = NULL, y = "Accuracy", color = "Method")
ggsave(file.path(out_dir, "boundary_interior_accuracy.png"), p_boundary, width = 12, height = 6, dpi = 220)

p_runtime <- summary_metrics |>
  filter(seconds_mean > 0) |>
  ggplot(aes(pattern, seconds_mean, color = method, group = method)) +
  geom_point(position = position_dodge(width = 0.35), size = 2) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = NULL, y = "Runtime, seconds (log scale)", color = "Method")
ggsave(file.path(out_dir, "geometric_runtime.png"), p_runtime, width = 10, height = 5, dpi = 220)

map_data <- lapply(seq_along(patterns), function(i) {
  sim <- simulate_spatial_domains(12000, patterns[i], k = 5, noise = 0.25, seed = 700 + i)
  refined <- refine_spatial_clusters(sim$xy, sim$labels)
  keep <- sample.int(nrow(sim$xy), min(6500L, nrow(sim$xy)))
  data.frame(
    x = sim$xy[keep, 1L], y = sim$xy[keep, 2L], pattern = patterns[i],
    Truth = sim$truth[keep], Initial = sim$labels[keep], Refined = refined[keep]
  ) |>
    pivot_longer(c(Truth, Initial, Refined), names_to = "state", values_to = "label")
}) |>
  bind_rows() |>
  mutate(state = factor(state, levels = c("Truth", "Initial", "Refined")))

p_atlas <- ggplot(map_data, aes(x, y, color = label)) +
  geom_point(size = 0.18, alpha = 0.8) +
  coord_equal() +
  facet_grid(pattern ~ state) +
  guides(color = "none") +
  theme(
    axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), strip.text.y = element_text(angle = 0)
  )
ggsave(file.path(out_dir, "geometric_shape_atlas.png"), p_atlas, width = 10, height = 18, dpi = 220)

scale_design <- tidyr::crossing(
  n = c(5000L, 20000L, 60000L, if (quick) integer() else c(150000L, 500000L)),
  dimensions = c(2L, 3L),
  replicate = seq_len(if (quick) 2L else 10L)
)
scale_results <- vector("list", nrow(scale_design))
for (i in seq_len(nrow(scale_design))) {
  row <- scale_design[i, ]
  pattern <- if (row$dimensions == 3L) "layers3d" else "jagged_stripes"
  sim <- simulate_spatial_domains(row$n, pattern, dimensions = row$dimensions, noise = 0.25, seed = 90000 + i)
  timing <- system.time(refined <- refine_spatial_clusters(sim$xy, sim$labels))
  scale_results[[i]] <- mutate(row, seconds = unname(timing["elapsed"]), accuracy = mean(refined == sim$truth))
}
scale_metrics <- bind_rows(scale_results)
write.csv(scale_metrics, file.path(out_dir, "scaling_metrics.csv"), row.names = FALSE)

p_scale <- ggplot(scale_metrics, aes(n, seconds, color = factor(dimensions))) +
  geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Observations", y = "Runtime, seconds", color = "Dimensions")
ggsave(file.path(out_dir, "runtime_scaling_graph.png"), p_scale, width = 7, height = 5, dpi = 220)

message("Wrote publication benchmark outputs to ", out_dir)
