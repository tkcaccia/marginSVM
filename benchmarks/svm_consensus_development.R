#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialGraphRefine)
  library(dplyr)
  library(ggplot2)
})

out_dir <- file.path("benchmarks", "results", "svm_upgrade")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workers <- max(1L, as.integer(Sys.getenv("SPATIAL_REFINE_BENCHMARK_WORKERS", "4")))

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
thresholds <- c(0, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20)
rows <- list()
position <- 0L
for (scenario in seq_len(length(cases) * 3L)) {
  case <- rep(cases, each = 3L)[scenario]
  replicate <- rep(1:3, times = length(cases))[scenario]
  sim <- make_case(case, 350000L + scenario)
  pred <- refine_spatial_svm(
    sim$xy, sim$labels, sim$samples, tiles = "auto", seed = 360000L + scenario,
    execution = list(
      target_tile_size = 4000L, overlap = 0.50, workers = workers,
      blend = "soft", adaptive = "multiscale"
    )
  )
  for (threshold in thresholds) {
    gated <- pred
    retain <- attr(pred, "margin") < threshold
    gated[retain] <- sim$labels[retain]
    position <- position + 1L
    rows[[position]] <- data.frame(
      case = case, replicate = replicate, threshold = threshold,
      accuracy = mean(gated == sim$truth), changed = mean(gated != sim$labels),
      correction = mean(gated[sim$labels != sim$truth] == sim$truth[sim$labels != sim$truth]),
      damage = mean(gated[sim$labels == sim$truth] != sim$truth[sim$labels == sim$truth])
    )
  }
}
metrics <- bind_rows(rows)
write.csv(metrics, file.path(out_dir, "consensus_development_metrics.csv"), row.names = FALSE)
summary <- metrics |>
  group_by(threshold) |>
  summarise(
    accuracy = mean(accuracy), worst_case = min(accuracy),
    correction = mean(correction), damage = mean(damage), changed = mean(changed),
    .groups = "drop"
  ) |>
  arrange(desc(accuracy), damage)
write.csv(summary, file.path(out_dir, "consensus_development_selection.csv"), row.names = FALSE)

p <- metrics |>
  group_by(case, threshold) |>
  summarise(accuracy = mean(accuracy), .groups = "drop") |>
  ggplot(aes(threshold, accuracy, color = case)) +
  geom_line() + geom_point() +
  labs(x = "Minimum multiscale margin required to change a label",
       y = "Development accuracy", color = "Structure") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "consensus_margin_development.png"), p,
       width = 9, height = 5.5, dpi = 220, bg = "white")

message("Consensus-gate development complete: ", out_dir)
