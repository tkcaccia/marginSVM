#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialGraphRefine)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mclust)
})

out_dir <- file.path("benchmarks", "results", "svm_gradient")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workers <- max(1L, as.integer(Sys.getenv("SPATIAL_REFINE_BENCHMARK_WORKERS", "4")))
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank()))

safe_mean <- function(x) if (length(x)) mean(x) else NA_real_

score <- function(pred, sim, method, elapsed) {
  correct <- pred == sim$truth
  initially_wrong <- sim$labels != sim$truth
  initially_correct <- !initially_wrong
  recalls <- vapply(levels(sim$truth), function(level) {
    safe_mean(pred[sim$truth == level] == level)
  }, numeric(1L))
  data.frame(
    method = method,
    accuracy = mean(correct),
    ari = mclust::adjustedRandIndex(pred, sim$truth),
    macro_recall = mean(recalls),
    worst_recall = min(recalls),
    boundary_accuracy = safe_mean(correct[sim$boundary_proximity <= 0.03]),
    interior_accuracy = safe_mean(correct[sim$boundary_proximity > 0.03]),
    correction_recall = safe_mean(correct[initially_wrong]),
    damage_rate = safe_mean(!correct[initially_correct]),
    seconds = unname(elapsed),
    stringsAsFactors = FALSE
  )
}

svm_grid <- crossing(
  gamma = c(8, 16, 32),
  n_features = c(128L, 256L),
  epochs = c(8L, 16L),
  learning_rate = c(0.005, 0.01)
)
development_design <- crossing(
  minority = c(0.05, 0.25, 0.40),
  replicate = 1:3
)

message("Selecting SVM settings on development seeds")
run_development <- function(index) {
  scenario <- development_design[index, ]
  sim <- simulate_gradient_regions(
    n = 12000, minority = scenario$minority, dimensions = 2,
    samples = 1, seed = 110000L + index
  )
  bind_rows(lapply(seq_len(nrow(svm_grid)), function(candidate) {
    settings <- svm_grid[candidate, ]
    timing <- system.time(pred <- refine_spatial_svm(
      sim$xy, sim$labels, sim$samples, tiles = NULL,
      gamma = settings$gamma, n_features = settings$n_features,
      epochs = settings$epochs, learning_rate = settings$learning_rate,
      seed = 120000L + candidate
    ))
    data.frame(scenario, candidate, settings,
               accuracy = mean(pred == sim$truth), seconds = unname(timing["elapsed"]))
  }))
}
development <- if (.Platform$OS.type != "windows" && workers > 1L) {
  parallel::mclapply(seq_len(nrow(development_design)), run_development,
                     mc.cores = workers, mc.preschedule = TRUE)
} else {
  lapply(seq_len(nrow(development_design)), run_development)
}
development_metrics <- bind_rows(development)
development_summary <- development_metrics |>
  group_by(candidate, gamma, n_features, epochs, learning_rate) |>
  summarise(accuracy = mean(accuracy), seconds = mean(seconds), .groups = "drop") |>
  arrange(desc(accuracy), seconds)
selected <- development_summary[1L, ]
write.csv(development_metrics, file.path(out_dir, "development_metrics.csv"), row.names = FALSE)
write.csv(development_summary, file.path(out_dir, "development_selection.csv"), row.names = FALSE)
write.csv(selected, file.path(out_dir, "selected_svm_settings.csv"), row.names = FALSE)
message(
  "Selected gamma=", selected$gamma, ", features=", selected$n_features,
  ", epochs=", selected$epochs, ", learning_rate=", selected$learning_rate
)

svm_prediction <- function(sim) {
  refine_spatial_svm(
    sim$xy, sim$labels, sim$samples, tiles = "auto",
    gamma = selected$gamma, n_features = selected$n_features,
    epochs = selected$epochs, learning_rate = selected$learning_rate,
    seed = 130001L,
    execution = list(target_tile_size = 4000L, overlap = 0.50, workers = 1L)
  )
}

methods <- list(
  "SVM refinement" = svm_prediction,
  "Graph refinement" = function(sim) refine_spatial_clusters(sim$xy, sim$labels, sim$samples),
  "SpaGCN refine" = function(sim) SpatialGraphRefine:::.refine_published_labels(
    sim$xy, sim$labels, sim$samples, method = "spagcn", neighbors = 6L
  ),
  "GraphST refine" = function(sim) SpatialGraphRefine:::.refine_published_labels(
    sim$xy, sim$labels, sim$samples, method = "graphst", neighbors = 50L
  ),
  "C++ kNN vote" = function(sim) refine_spatial_clusters(
    sim$xy, sim$labels, sim$samples,
    control = list(weighted = FALSE, iterations = 1L, consensus = 0.5,
                   preserve = 0, margin = 0, current_support = 1)
  ),
  "Potts-like ICM" = function(sim) refine_spatial_clusters(
    sim$xy, sim$labels, sim$samples,
    control = list(weighted = FALSE, iterations = 8L, consensus = 0.5,
                   preserve = 0.35, margin = 0, current_support = 1)
  )
)

test_design <- crossing(
  minority = c(0.05, 0.15, 0.25, 0.35, 0.40),
  dimensions = c(2L, 3L),
  tissues = c(1L, 3L),
  replicate = 4:10
)

message("Running held-out gradient-region benchmark")
run_test <- function(index) {
  scenario <- test_design[index, ]
  sim <- simulate_gradient_regions(
    n = 12000, minority = scenario$minority,
    dimensions = scenario$dimensions, samples = scenario$tissues,
    seed = 140000L + index
  )
  rows <- list(score(sim$labels, sim, "Initial", 0))
  for (method in names(methods)) {
    timing <- system.time(pred <- methods[[method]](sim))
    rows[[length(rows) + 1L]] <- score(pred, sim, method, timing["elapsed"])
  }
  bind_cols(scenario[rep(1L, length(rows)), ], bind_rows(rows))
}
test <- if (.Platform$OS.type != "windows" && workers > 1L) {
  parallel::mclapply(seq_len(nrow(test_design)), run_test,
                     mc.cores = workers, mc.preschedule = TRUE)
} else {
  lapply(seq_len(nrow(test_design)), run_test)
}
test_metrics <- bind_rows(test)
write.csv(test_metrics, file.path(out_dir, "heldout_metrics.csv"), row.names = FALSE)

summary_metrics <- test_metrics |>
  group_by(method) |>
  summarise(
    accuracy = mean(accuracy), ari = mean(ari), macro_recall = mean(macro_recall),
    worst_recall = mean(worst_recall), boundary_accuracy = mean(boundary_accuracy),
    interior_accuracy = mean(interior_accuracy), correction_recall = mean(correction_recall),
    damage_rate = mean(damage_rate), seconds = mean(seconds), .groups = "drop"
  ) |>
  arrange(desc(accuracy))
write.csv(summary_metrics, file.path(out_dir, "heldout_summary.csv"), row.names = FALSE)

paired <- test_metrics |>
  select(minority, dimensions, tissues, replicate, method, accuracy) |>
  pivot_wider(names_from = method, values_from = accuracy)
set.seed(150001L)
comparators <- setdiff(names(methods), "SVM refinement")
paired_intervals <- bind_rows(lapply(comparators, function(comparator) {
  difference <- paired[["SVM refinement"]] - paired[[comparator]]
  bootstrap <- replicate(10000L, mean(sample(difference, replace = TRUE)))
  data.frame(
    comparator = comparator, mean_difference = mean(difference),
    lower_95 = unname(quantile(bootstrap, 0.025)),
    upper_95 = unname(quantile(bootstrap, 0.975))
  )
}))
write.csv(paired_intervals, file.path(out_dir, "svm_paired_effects.csv"), row.names = FALSE)

by_condition <- test_metrics |>
  group_by(minority, dimensions, tissues, method) |>
  summarise(accuracy = mean(accuracy), damage_rate = mean(damage_rate), .groups = "drop")
p_accuracy <- ggplot(by_condition, aes(minority, accuracy, color = method)) +
  geom_line() + geom_point(size = 1.2) +
  facet_grid(dimensions ~ tissues, labeller = label_both) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(x = "Minority-label fraction per area", y = "Held-out reference accuracy", color = "Method") +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "gradient_accuracy.png"), p_accuracy,
       width = 11, height = 6.5, dpi = 220)

advantages <- by_condition |>
  select(minority, dimensions, tissues, method, accuracy) |>
  pivot_wider(names_from = method, values_from = accuracy) |>
  mutate(
    svm_minus_graph = `SVM refinement` - `Graph refinement`,
    svm_minus_spagcn = `SVM refinement` - `SpaGCN refine`,
    svm_minus_graphst = `SVM refinement` - `GraphST refine`
  )
write.csv(advantages, file.path(out_dir, "svm_advantage_by_condition.csv"), row.names = FALSE)

message("Running 60k multi-tissue overlapping-tile example")
large <- simulate_gradient_regions(
  n = 60000, minority = 0.25, dimensions = 3, samples = 3, seed = 160001L
)
untiled_time <- system.time(untiled <- refine_spatial_svm(
  large$xy, large$labels, large$samples, tiles = NULL,
  gamma = selected$gamma, n_features = selected$n_features,
  epochs = selected$epochs, learning_rate = selected$learning_rate, seed = 160002L
))
tiled_time <- system.time(tiled <- refine_spatial_svm(
  large$xy, large$labels, large$samples, tiles = "auto",
  gamma = selected$gamma, n_features = selected$n_features,
  epochs = selected$epochs, learning_rate = selected$learning_rate, seed = 160002L,
  execution = list(target_tile_size = 10000L, overlap = 0.50, workers = workers)
))
large_summary <- bind_rows(
  data.frame(
    method = "Untiled SVM", accuracy = mean(untiled == large$truth),
    seconds = unname(untiled_time["elapsed"]), workers = 1L,
    tile_agreement = 1
  ),
  data.frame(
    method = "Overlapping tiled SVM", accuracy = mean(tiled == large$truth),
    seconds = unname(tiled_time["elapsed"]), workers = attr(tiled, "workers"),
    tile_agreement = mean(tiled == untiled)
  )
)
write.csv(large_summary, file.path(out_dir, "multitissue_60k_scaling.csv"), row.names = FALSE)
write.csv(attr(tiled, "tiles"), file.path(out_dir, "multitissue_60k_tiles.csv"), row.names = FALSE)

message("Generating requested four-area example")
example <- simulate_gradient_regions(n = 50000, minority = 0.05, seed = 170001L)
example_svm <- svm_prediction(example)
example_graph <- refine_spatial_clusters(example$xy, example$labels, example$samples)
area_purity <- bind_rows(lapply(
  list("Mixed input" = example$labels, "SVM refinement" = example_svm,
       "Graph refinement" = example_graph),
  function(pred) bind_rows(lapply(levels(example$area), function(area) {
    rows <- example$area == area
    reference <- unique(example$truth[rows])
    data.frame(area = area, reference = as.character(reference), purity = mean(pred[rows] == reference))
  }))
), .id = "method")
write.csv(area_purity, file.path(out_dir, "requested_example_area_purity.csv"), row.names = FALSE)

set.seed(170002L)
display <- sample.int(nrow(example$xy), 35000L)
map_data <- bind_rows(
  data.frame(example$xy[display, ], panel = "Reference A-B-B-C", label = example$truth[display]),
  data.frame(example$xy[display, ], panel = "95%-5% mixed labels", label = example$labels[display]),
  data.frame(example$xy[display, ], panel = "SVM refinement", label = example_svm[display]),
  data.frame(example$xy[display, ], panel = "Graph refinement", label = example_graph[display])
)
map_data$panel <- factor(
  map_data$panel,
  c("Reference A-B-B-C", "95%-5% mixed labels", "SVM refinement", "Graph refinement")
)
p_map <- ggplot(map_data, aes(x, y, color = label)) +
  geom_point(size = 0.25, alpha = 0.8) + facet_wrap(~panel, ncol = 2) +
  coord_equal() + labs(x = NULL, y = NULL, color = "Group") +
  theme_void(base_size = 10) +
  theme(
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
ggsave(file.path(out_dir, "requested_gradient_example.png"), p_map,
       width = 8.5, height = 7, dpi = 240, bg = "white")

message("SVM gradient benchmark complete: ", out_dir)
