#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(marginSVM)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

out_dir <- file.path("benchmarks", "results", "marginsvm_speed_tuning")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workers <- max(1L, as.integer(Sys.getenv("SPATIAL_REFINE_BENCHMARK_WORKERS", "4")))

design <- crossing(
  pattern = c("jagged_stripes", "rings", "thin_layers", "islands", "intermixed", "layers3d"),
  density_profile = c("uniform", "extreme"),
  noise_type = c("random", "boundary"),
  replicate = 1:2
) |>
  mutate(scenario = row_number(), dimensions = if_else(pattern == "layers3d", 3L, 2L))

settings <- list(
  "Current default" = list(),
  "Halo 0.25" = list(overlap = 0.25),
  "Halo 0.15" = list(overlap = 0.15),
  "48 landmarks" = list(landmarks = 48L),
  "8 epochs" = list(epochs = 8L),
  "24 TV iterations" = list(tv_iterations = 24L),
  "Fast balanced" = list(overlap = 0.25, landmarks = 48L, epochs = 8L, tv_iterations = 24L),
  "Conservative fast" = list(overlap = 0.25, landmarks = 64L, epochs = 8L, tv_iterations = 24L)
)

run_one <- function(index) {
  d <- design[index, ]
  sim <- simulate_spatial_domains(
    12000L, d$pattern, dimensions = d$dimensions, noise = 0.25,
    noise_type = d$noise_type, density_profile = d$density_profile,
    seed = 3100000L + index
  )
  initially_wrong <- sim$labels != sim$truth
  boundary <- sim$boundary_proximity <= quantile(sim$boundary_proximity, 0.2)
  result <- bind_rows(lapply(seq_along(settings), function(j) {
    control <- utils::modifyList(
      list(workers = 1L, seed = 3200000L + index), settings[[j]])
    timing <- system.time(pred <- marginSVM:::.refine_spatial_svm_engine(
      sim$xy, sim$labels, sim$samples, control = control))
    correct <- pred == sim$truth
    data.frame(
      setting = names(settings)[j],
      accuracy = mean(correct),
      boundary_accuracy = mean(correct[boundary]),
      correction_recall = mean(correct[initially_wrong]),
      damage_rate = mean(!correct[!initially_wrong]),
      seconds = unname(timing[["elapsed"]])
    )
  }))
  bind_cols(d[rep(1L, length(settings)), ], result)
}

message("Running ", nrow(design), " paired speed-quality scenarios")
rows <- if (.Platform$OS.type != "windows" && workers > 1L) {
  parallel::mclapply(seq_len(nrow(design)), run_one, mc.cores = workers,
                     mc.preschedule = TRUE)
} else {
  lapply(seq_len(nrow(design)), run_one)
}
metrics <- bind_rows(rows)
write.csv(metrics, file.path(out_dir, "speed_tuning_metrics.csv"), row.names = FALSE)

baseline <- metrics |>
  filter(setting == "Current default") |>
  select(scenario, baseline_accuracy = accuracy, baseline_seconds = seconds)
paired <- metrics |>
  left_join(baseline, by = "scenario") |>
  group_by(setting) |>
  summarise(
    accuracy = mean(accuracy),
    accuracy_delta = mean(accuracy - baseline_accuracy),
    boundary_accuracy = mean(boundary_accuracy),
    correction_recall = mean(correction_recall),
    damage_rate = mean(damage_rate),
    seconds = mean(seconds),
    speedup = mean(baseline_seconds) / mean(seconds),
    .groups = "drop"
  ) |>
  arrange(desc(speedup))
write.csv(paired, file.path(out_dir, "speed_tuning_summary.csv"), row.names = FALSE)

plot_data <- paired |>
  filter(setting != "Current default") |>
  mutate(eligible = accuracy_delta >= -0.001)
plot <- ggplot(plot_data, aes(speedup, accuracy_delta, label = setting, color = eligible)) +
  geom_hline(yintercept = -0.001, color = "grey55", linetype = 2, linewidth = 0.5) +
  geom_vline(xintercept = 1, color = "grey75", linewidth = 0.5) +
  geom_point(size = 2.6) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf, seed = 1) +
  scale_color_manual(values = c(`TRUE` = "#2A9D8F", `FALSE` = "#C44E52"), guide = "none") +
  labs(x = "Speedup relative to current default", y = "Paired accuracy difference") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "speed_accuracy_tradeoff.png"), plot,
       width = 8.2, height = 5.0, dpi = 240, bg = "white")

print(paired)
