#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialGraphRefine)
  library(ggplot2)
  library(mclust)
})

out_dir <- file.path("benchmarks", "results", "structured_svm_geometry")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
patterns <- c(
  "jagged_stripes", "wavy_layers", "rings", "spiral", "branching",
  "lobes", "islands", "disconnected", "thin_layers", "intermixed",
  "layers3d"
)
replicates <- 1:3
methods <- c("Structured SVM", "Graph refinement", "GraphST refine",
             "SpaGCN refine", "Potts-like ICM")

run_methods <- function(sim) {
  list(
    "Structured SVM" = function() refine_spatial_svm(sim$xy, sim$labels),
    "Graph refinement" = function() refine_spatial_clusters(sim$xy, sim$labels),
    "GraphST refine" = function() SpatialGraphRefine:::.refine_published_labels(
      sim$xy, sim$labels, sim$samples, method = "graphst", neighbors = 50L
    ),
    "SpaGCN refine" = function() SpatialGraphRefine:::.refine_published_labels(
      sim$xy, sim$labels, sim$samples, method = "spagcn", neighbors = 6L
    ),
    "Potts-like ICM" = function() refine_spatial_clusters(
      sim$xy, sim$labels, sim$samples,
      control = list(weighted = FALSE, iterations = 8L, consensus = 0.5,
                     preserve = 0.35, margin = 0, current_support = 1)
    )
  )
}

records <- list()
example <- NULL
index <- 0L
for (pattern in patterns) {
  for (replicate in replicates) {
    sim <- simulate_spatial_domains(
      n = 10000L, pattern = pattern, k = 5L, noise = 0.25,
      noise_type = if (replicate == 1L) "random" else if (replicate == 2L) "boundary" else "region",
      seed = 730000L + match(pattern, patterns) * 100L + replicate
    )
    predictions <- list(Truth = sim$truth, Initial = sim$labels)
    boundary_rows <- sim$boundary_proximity <=
      stats::quantile(sim$boundary_proximity, 0.20, names = FALSE)
    for (method in names(run_methods(sim))) {
      timing <- system.time(prediction <- run_methods(sim)[[method]]())
      predictions[[method]] <- prediction
      index <- index + 1L
      records[[index]] <- data.frame(
        pattern = pattern, replicate = replicate,
        noise_type = c("random", "boundary", "region")[replicate],
        method = method,
        accuracy = mean(prediction == sim$truth),
        ari = adjustedRandIndex(prediction, sim$truth),
        boundary_accuracy = mean(prediction[boundary_rows] == sim$truth[boundary_rows]),
        seconds = unname(timing[["elapsed"]])
      )
    }
    if (pattern == "jagged_stripes" && replicate == 2L) {
      example <- list(sim = sim, predictions = predictions)
    }
  }
}

results <- do.call(rbind, records)
results$method <- factor(results$method, levels = methods)
write.csv(results, file.path(out_dir, "heldout_metrics.csv"), row.names = FALSE)
summary <- aggregate(cbind(accuracy, ari, boundary_accuracy, seconds) ~ method,
                     results, mean)
write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)

p_accuracy <- ggplot(results, aes(method, accuracy, colour = method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12), size = 0.9, alpha = 0.65) +
  facet_wrap(~pattern, ncol = 5) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none") +
  labs(x = NULL, y = "Assignment accuracy")
ggsave(file.path(out_dir, "accuracy_by_geometry.png"), p_accuracy,
       width = 13, height = 7, dpi = 220)

plot_rows <- lapply(names(example$predictions), function(name) {
  data.frame(
    x = example$sim$xy[, 1L], y = example$sim$xy[, 2L],
    label = example$predictions[[name]], panel = name
  )
})
plot_data <- do.call(rbind, plot_rows)
plot_data$panel <- factor(plot_data$panel, levels = names(example$predictions))
p_example <- ggplot(plot_data, aes(x, y, colour = label)) +
  geom_point(size = 0.18, alpha = 0.9) +
  coord_equal() +
  facet_wrap(~panel, ncol = 4) +
  theme_void(base_size = 10) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave(file.path(out_dir, "jagged_stripes_predictions.png"), p_example,
       width = 12, height = 6.8, dpi = 240)

print(summary[order(-summary$accuracy), ])
