#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialGraphRefine)
  library(ggplot2)
})

out_dir <- file.path("benchmarks", "results", "structured_svm_gradient")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sim <- simulate_gradient_regions(
  n = 50000L, minority = 0.05, samples = 2L, seed = 7L
)

elapsed <- system.time(
  structured <- refine_spatial_svm(sim$xy, sim$labels, sim$samples)
)[["elapsed"]]
graph <- refine_spatial_clusters(sim$xy, sim$labels, sim$samples)

area_accuracy <- vapply(levels(sim$area), function(area) {
  rows <- sim$area == area
  mean(structured[rows] == sim$truth[rows])
}, numeric(1L))
metrics <- data.frame(
  method = c("Initial", "Structured SVM", "Graph refinement"),
  accuracy = c(mean(sim$labels == sim$truth), mean(structured == sim$truth),
               mean(graph == sim$truth)),
  seconds = c(0, elapsed, NA_real_)
)
write.csv(metrics, file.path(out_dir, "metrics.csv"), row.names = FALSE)
write.csv(data.frame(area = names(area_accuracy), accuracy = area_accuracy),
          file.path(out_dir, "area_accuracy.csv"), row.names = FALSE)

set.seed(8L)
display <- sample.int(nrow(sim$xy), 30000L)
panels <- list(
  Reference = sim$truth,
  `Mixed input` = sim$labels,
  `Structured SVM` = structured,
  `Graph refinement` = graph
)
plot_data <- do.call(rbind, lapply(names(panels), function(panel) {
  data.frame(
    x = sim$xy[display, 1L], y = sim$xy[display, 2L],
    tissue = sim$samples[display], label = panels[[panel]][display], panel = panel
  )
}))
plot_data$panel <- factor(plot_data$panel, levels = names(panels))
p <- ggplot(plot_data, aes(x, y, colour = label)) +
  geom_point(size = 0.18, alpha = 0.85) +
  coord_equal() +
  facet_wrap(~panel, ncol = 2) +
  theme_void(base_size = 10) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave(file.path(out_dir, "gradient_50k.png"), p,
       width = 9.5, height = 7.5, dpi = 240)

print(metrics)
print(area_accuracy)
