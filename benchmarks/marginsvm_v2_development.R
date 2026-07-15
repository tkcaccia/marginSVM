#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialGraphRefine)
  library(ggplot2)
  library(mclust)
})

out_dir <- file.path("benchmarks", "results", "marginsvm_v2_development")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workers <- max(1L, as.integer(Sys.getenv("SPATIAL_REFINE_BENCHMARK_WORKERS", "4")))
n <- as.integer(Sys.getenv("MARGINSVM_V2_N", "6000"))
replicates <- as.integer(Sys.getenv("MARGINSVM_V2_REPLICATES", "2"))

design <- expand.grid(
  pattern = c("jagged_stripes", "wavy_layers", "rings", "spiral",
              "branching", "islands", "thin_layers", "layers3d"),
  noise_type = c("random", "boundary"),
  replicate = seq_len(replicates),
  stringsAsFactors = FALSE
)
design <- subset(design,
  (noise_type == "boundary" & pattern %in% c("jagged_stripes", "wavy_layers", "spiral", "thin_layers")) |
  (noise_type == "random" & pattern %in% c("rings", "branching", "islands", "layers3d")))

variants <- list(
  "marginSVM v1" = list(experimental_v2 = 0),
  "marginSVM v2" = list(experimental_v2 = 1),
  "v2 isotropic" = list(experimental_v2 = 1, anisotropy = 0),
  "v2 no specialists" = list(experimental_v2 = 1, pairwise_specialists = 0)
)

rows <- list()
examples <- list()
position <- 0L
example_position <- 0L
for (scenario in seq_len(nrow(design))) {
  d <- design[scenario, ]
  sim <- simulate_spatial_domains(
    n = n, pattern = d$pattern, dimensions = if (d$pattern == "layers3d") 3L else 2L,
    k = 5L, noise = 0.25, noise_type = d$noise_type,
    seed = 900000L + scenario
  )
  boundary <- sim$boundary_proximity <= 0.08
  for (method in names(variants)) {
    control <- c(variants[[method]], list(workers = workers, seed = 910000L + scenario))
    timing <- system.time(pred <- refine_spatial_svm(
      sim$xy, sim$labels, sim$samples, control = control
    ))
    decision <- attr(pred, "decision")
    position <- position + 1L
    rows[[position]] <- data.frame(
      scenario = scenario, pattern = d$pattern, noise_type = d$noise_type,
      replicate = d$replicate, method = method,
      accuracy = mean(pred == sim$truth),
      ari = adjustedRandIndex(pred, sim$truth),
      boundary_accuracy = mean(pred[boundary] == sim$truth[boundary]),
      rare_accuracy = mean(pred[sim$truth != levels(sim$truth)[1L]] ==
                           sim$truth[sim$truth != levels(sim$truth)[1L]]),
      damage = mean(pred[!sim$corrupted] != sim$truth[!sim$corrupted]),
      corrected = mean(pred[sim$corrupted] == sim$truth[sim$corrupted]),
      changed = if (is.null(decision)) mean(pred != sim$labels) else mean(decision == "change"),
      unresolved = if (is.null(decision)) 0 else mean(decision == "unresolved"),
      seconds = unname(timing["elapsed"]), stringsAsFactors = FALSE
    )
    if (d$replicate == 1L && d$pattern %in% c("jagged_stripes", "thin_layers") &&
        method %in% c("marginSVM v1", "marginSVM v2")) {
      display <- if (n > 12000L) sample.int(n, 12000L) else seq_len(n)
      example_position <- example_position + 1L
      examples[[example_position]] <- data.frame(
        x = sim$xy[display, 1L], y = sim$xy[display, 2L], label = pred[display],
        pattern = d$pattern, panel = method, stringsAsFactors = FALSE
      )
      if (method == "marginSVM v1") {
        example_position <- example_position + 1L
        examples[[example_position]] <- data.frame(
          x = sim$xy[display, 1L], y = sim$xy[display, 2L], label = sim$truth[display],
          pattern = d$pattern, panel = "Reference", stringsAsFactors = FALSE
        )
        example_position <- example_position + 1L
        examples[[example_position]] <- data.frame(
          x = sim$xy[display, 1L], y = sim$xy[display, 2L], label = sim$labels[display],
          pattern = d$pattern, panel = "Corrupted input", stringsAsFactors = FALSE
        )
      }
    }
  }
  message("completed scenario ", scenario, "/", nrow(design))
}

metrics <- do.call(rbind, rows)
write.csv(metrics, file.path(out_dir, "development_metrics.csv"), row.names = FALSE)
summary <- aggregate(cbind(accuracy, ari, boundary_accuracy, rare_accuracy, damage,
                           corrected, changed, unresolved, seconds) ~ method,
                     metrics, mean)
write.csv(summary, file.path(out_dir, "development_summary.csv"), row.names = FALSE)

long <- reshape(metrics,
  varying = c("accuracy", "boundary_accuracy", "rare_accuracy"),
  v.names = "value", timevar = "metric",
  times = c("Overall", "Boundary band", "Non-background / rare"), direction = "long")
long$method <- factor(long$method, names(variants))
p_metrics <- ggplot(long, aes(method, value, fill = method)) +
  geom_boxplot(width = 0.68, outlier.size = 0.8) +
  facet_wrap(~metric, nrow = 1L) +
  coord_cartesian(ylim = c(0.45, 1)) +
  scale_fill_manual(values = c("#4C78A8", "#E45756", "#72B7B2", "#F2CF5B")) +
  labs(x = NULL, y = "Reference accuracy") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none",
        panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "development_accuracy.png"), p_metrics,
       width = 10.5, height = 4.2, dpi = 240, bg = "white")

map_data <- do.call(rbind, examples)
map_data$panel <- factor(map_data$panel,
  c("Reference", "Corrupted input", "marginSVM v1", "marginSVM v2"))
p_maps <- ggplot(map_data, aes(x, y, color = label)) +
  geom_point(size = 0.22, alpha = 0.85) + coord_equal() +
  facet_grid(pattern ~ panel) +
  guides(color = guide_legend(nrow = 1L, override.aes = list(size = 2))) +
  labs(x = NULL, y = NULL, color = "Region") + theme_void(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
        plot.background = element_rect(fill = "white", color = NA))
ggsave(file.path(out_dir, "development_spatial_examples.png"), p_maps,
       width = 12, height = 6.2, dpi = 240, bg = "white")

print(summary)
