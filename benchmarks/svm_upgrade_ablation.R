#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialGraphRefine)
  library(dplyr)
  library(ggplot2)
  library(mclust)
})

out_dir <- file.path("benchmarks", "results", "svm_upgrade")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workers <- max(1L, as.integer(Sys.getenv("SPATIAL_REFINE_BENCHMARK_WORKERS", "4")))
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank()))

variants <- list(
  "Hard regular tiles" = list(blend = "hard", adaptive = FALSE),
  "Soft regular tiles" = list(blend = "soft", adaptive = FALSE),
  "Adaptive soft tiles" = list(blend = "soft", adaptive = TRUE),
  "Multiscale soft ensemble" = list(blend = "soft", adaptive = "multiscale")
)

make_case <- function(case, seed, n = 20000L) {
  if (grepl("gradient", case, fixed = TRUE)) {
    minority <- if (grepl("40", case, fixed = TRUE)) 0.40 else 0.25
    return(simulate_gradient_regions(n, minority = minority, seed = seed))
  }
  pieces <- strsplit(case, "_")[[1L]]
  noise <- as.numeric(tail(pieces, 1L)) / 100
  pattern <- paste(head(pieces, -1L), collapse = "_")
  simulate_spatial_domains(n, pattern = pattern, k = 5L, noise = noise, seed = seed)
}

cases <- c("gradient_25", "gradient_40", "jagged_stripes_25", "rings_25", "lobes_25")
design <- expand.grid(case = cases, replicate = 1:3, stringsAsFactors = FALSE)

rows <- list()
position <- 0L
for (scenario in seq_len(nrow(design))) {
  sim <- make_case(design$case[scenario], 310000L + scenario)
  for (variant in names(variants)) {
    settings <- variants[[variant]]
    timing <- system.time(pred <- refine_spatial_svm(
      sim$xy, sim$labels, sim$samples, tiles = "auto", seed = 320000L + scenario,
      execution = list(
        target_tile_size = 4000L, overlap = 0.50, workers = workers,
        blend = settings$blend, adaptive = settings$adaptive
      )
    ))
    position <- position + 1L
    rows[[position]] <- data.frame(
      case = design$case[scenario], replicate = design$replicate[scenario],
      method = variant, accuracy = mean(pred == sim$truth),
      ari = mclust::adjustedRandIndex(pred, sim$truth),
      boundary_accuracy = mean(pred[sim$boundary_proximity <= 0.04] ==
                                 sim$truth[sim$boundary_proximity <= 0.04]),
      seconds = unname(timing["elapsed"]), tiles = nrow(attr(pred, "tiles")),
      mean_margin = mean(attr(pred, "margin")), stringsAsFactors = FALSE
    )
  }
}
metrics <- bind_rows(rows)
write.csv(metrics, file.path(out_dir, "upgrade_ablation_metrics.csv"), row.names = FALSE)
summary <- metrics |>
  group_by(case, method) |>
  summarise(
    accuracy = mean(accuracy), boundary_accuracy = mean(boundary_accuracy),
    ari = mean(ari), seconds = mean(seconds), tiles = mean(tiles), .groups = "drop"
  )
write.csv(summary, file.path(out_dir, "upgrade_ablation_summary.csv"), row.names = FALSE)

p_ablation <- summary |>
  tidyr::pivot_longer(c(accuracy, boundary_accuracy), names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric, accuracy = "All observations",
                         boundary_accuracy = "Boundary band")) |>
  ggplot(aes(method, value, fill = method)) +
  geom_col(width = 0.72) +
  facet_grid(metric ~ case) +
  coord_cartesian(ylim = c(0.45, 1)) +
  labs(x = NULL, y = "Reference accuracy", fill = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom")
ggsave(file.path(out_dir, "upgrade_ablation_accuracy.png"), p_ablation,
       width = 12, height = 6.5, dpi = 220, bg = "white")

examples <- c("gradient_40", "jagged_stripes_25", "rings_25")
maps <- list()
map_position <- 0L
for (index in seq_along(examples)) {
  sim <- make_case(examples[index], 330000L + index, n = 50000L)
  display <- sample.int(nrow(sim$xy), 30000L)
  map_position <- map_position + 1L
  maps[[map_position]] <- data.frame(
    sim$xy[display, 1:2, drop = FALSE], case = examples[index], method = "Reference",
    label = sim$truth[display]
  )
  map_position <- map_position + 1L
  maps[[map_position]] <- data.frame(
    sim$xy[display, 1:2, drop = FALSE], case = examples[index], method = "Noisy input",
    label = sim$labels[display]
  )
  for (variant in names(variants)) {
    settings <- variants[[variant]]
    pred <- refine_spatial_svm(
      sim$xy, sim$labels, sim$samples, tiles = "auto", seed = 340000L + index,
      execution = list(
        target_tile_size = 6000L, overlap = 0.50, workers = workers,
        blend = settings$blend, adaptive = settings$adaptive
      )
    )
    map_position <- map_position + 1L
    maps[[map_position]] <- data.frame(
      sim$xy[display, 1:2, drop = FALSE], case = examples[index], method = variant,
      label = pred[display]
    )
  }
}
map_data <- bind_rows(maps)
map_data$method <- factor(
  map_data$method,
  c("Reference", "Noisy input", "Hard regular tiles", "Soft regular tiles",
    "Adaptive soft tiles", "Multiscale soft ensemble")
)
p_maps <- ggplot(map_data, aes(x, y, color = label)) +
  geom_point(size = 0.12, alpha = 0.75) +
  facet_grid(case ~ method) + coord_equal() +
  labs(x = NULL, y = NULL, color = "Region") + theme_void(base_size = 9) +
  theme(
    legend.position = "bottom", plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
ggsave(file.path(out_dir, "upgrade_spatial_examples.png"), p_maps,
       width = 14, height = 8.5, dpi = 240, bg = "white")

message("SVM upgrade ablation complete: ", out_dir)
