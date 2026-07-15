#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialGraphRefine)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mclust)
})

out_dir <- file.path("benchmarks", "results", "svm_geometric_heldout")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workers <- max(1L, as.integer(Sys.getenv("SPATIAL_REFINE_BENCHMARK_WORKERS", "4")))
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank()))

methods <- list(
  "SVM refinement" = function(sim) refine_spatial_svm(
    sim$xy, sim$labels, sim$samples, tiles = "auto", seed = 410001L,
    execution = list(target_tile_size = 4000L, overlap = 0.50, workers = 1L)
  ),
  "Graph refinement" = function(sim) refine_spatial_clusters(sim$xy, sim$labels, sim$samples),
  "SpaGCN refine" = function(sim) SpatialGraphRefine:::.refine_published_labels(
    sim$xy, sim$labels, sim$samples, method = "spagcn", neighbors = 6L
  ),
  "GraphST refine" = function(sim) SpatialGraphRefine:::.refine_published_labels(
    sim$xy, sim$labels, sim$samples, method = "graphst", neighbors = 50L
  ),
  "Potts-like ICM" = function(sim) refine_spatial_clusters(
    sim$xy, sim$labels, sim$samples,
    control = list(weighted = FALSE, iterations = 8L, consensus = 0.5,
                   preserve = 0.35, margin = 0, current_support = 1)
  )
)

design <- crossing(
  pattern = c("jagged_stripes", "wavy_layers", "rings", "branching", "lobes"),
  noise = c(0.15, 0.25, 0.40), noise_type = c("random", "boundary"), replicate = 4:8
)

run_scenario <- function(index) {
  scenario <- design[index, ]
  sim <- simulate_spatial_domains(
    n = 12000L, pattern = scenario$pattern, k = 5L,
    noise = scenario$noise, noise_type = scenario$noise_type,
    seed = 420000L + index
  )
  initially_wrong <- sim$labels != sim$truth
  rows <- list(data.frame(
    method = "Initial", accuracy = mean(sim$labels == sim$truth),
    ari = adjustedRandIndex(sim$labels, sim$truth),
    boundary_accuracy = mean(sim$labels[sim$boundary_proximity <= 0.04] ==
                               sim$truth[sim$boundary_proximity <= 0.04]),
    correction_recall = 0, damage_rate = 0, seconds = 0
  ))
  for (method in names(methods)) {
    timing <- system.time(pred <- methods[[method]](sim))
    rows[[length(rows) + 1L]] <- data.frame(
      method = method, accuracy = mean(pred == sim$truth),
      ari = adjustedRandIndex(pred, sim$truth),
      boundary_accuracy = mean(pred[sim$boundary_proximity <= 0.04] ==
                                 sim$truth[sim$boundary_proximity <= 0.04]),
      correction_recall = mean(pred[initially_wrong] == sim$truth[initially_wrong]),
      damage_rate = mean(pred[!initially_wrong] != sim$truth[!initially_wrong]),
      seconds = unname(timing["elapsed"])
    )
  }
  bind_cols(scenario[rep(1L, length(rows)), ], bind_rows(rows))
}

message("Running ", nrow(design), " held-out geometric scenarios")
results <- if (.Platform$OS.type != "windows" && workers > 1L) {
  parallel::mclapply(seq_len(nrow(design)), run_scenario,
                     mc.cores = workers, mc.preschedule = TRUE)
} else {
  lapply(seq_len(nrow(design)), run_scenario)
}
metrics <- bind_rows(results)
write.csv(metrics, file.path(out_dir, "geometric_heldout_metrics.csv"), row.names = FALSE)
summary <- metrics |>
  group_by(method) |>
  summarise(
    accuracy = mean(accuracy), ari = mean(ari),
    boundary_accuracy = mean(boundary_accuracy),
    correction_recall = mean(correction_recall), damage_rate = mean(damage_rate),
    seconds = mean(seconds), .groups = "drop"
  ) |>
  arrange(desc(accuracy))
write.csv(summary, file.path(out_dir, "geometric_heldout_summary.csv"), row.names = FALSE)

paired <- metrics |>
  select(pattern, noise, noise_type, replicate, method, accuracy) |>
  pivot_wider(names_from = method, values_from = accuracy)
set.seed(430001L)
effects <- bind_rows(lapply(setdiff(names(methods), "SVM refinement"), function(comparator) {
  difference <- paired[["SVM refinement"]] - paired[[comparator]]
  bootstrap <- replicate(10000L, mean(sample(difference, replace = TRUE)))
  data.frame(
    comparator = comparator, difference = mean(difference),
    lower_95 = unname(quantile(bootstrap, 0.025)),
    upper_95 = unname(quantile(bootstrap, 0.975))
  )
}))
write.csv(effects, file.path(out_dir, "geometric_paired_effects.csv"), row.names = FALSE)

plot_data <- metrics |>
  group_by(pattern, noise, noise_type, method) |>
  summarise(accuracy = mean(accuracy), .groups = "drop")
p_accuracy <- ggplot(plot_data, aes(noise, accuracy, color = method)) +
  geom_line() + geom_point(size = 1) +
  facet_grid(noise_type ~ pattern) + coord_cartesian(ylim = c(0.45, 1)) +
  labs(x = "Injected label-error fraction", y = "Held-out reference accuracy", color = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "geometric_heldout_accuracy.png"), p_accuracy,
       width = 13, height = 6.5, dpi = 220, bg = "white")

examples <- c("jagged_stripes", "rings", "lobes")
map_rows <- list()
position <- 0L
for (index in seq_along(examples)) {
  sim <- simulate_spatial_domains(
    n = 50000L, pattern = examples[index], k = 5L,
    noise = 0.25, noise_type = "random", seed = 440000L + index
  )
  predictions <- list(
    "Reference" = sim$truth, "Noisy input" = sim$labels,
    "SVM refinement" = methods[["SVM refinement"]](sim),
    "Graph refinement" = methods[["Graph refinement"]](sim),
    "GraphST refine" = methods[["GraphST refine"]](sim)
  )
  set.seed(450000L + index)
  display <- sample.int(nrow(sim$xy), 32000L)
  for (method in names(predictions)) {
    position <- position + 1L
    map_rows[[position]] <- data.frame(
      sim$xy[display, 1:2, drop = FALSE], pattern = examples[index], method = method,
      label = predictions[[method]][display]
    )
  }
}
map_data <- bind_rows(map_rows)
map_data$method <- factor(
  map_data$method,
  c("Reference", "Noisy input", "SVM refinement", "Graph refinement", "GraphST refine")
)
p_maps <- ggplot(map_data, aes(x, y, color = label)) +
  geom_point(size = 0.12, alpha = 0.75) + facet_grid(pattern ~ method) +
  coord_equal() + labs(x = NULL, y = NULL, color = "Region") + theme_void(base_size = 9) +
  theme(
    legend.position = "bottom", plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
ggsave(file.path(out_dir, "geometric_spatial_examples.png"), p_maps,
       width = 13, height = 8.5, dpi = 240, bg = "white")

message("Held-out geometric benchmark complete: ", out_dir)
