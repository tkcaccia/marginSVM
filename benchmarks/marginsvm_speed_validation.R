#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialGraphRefine)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

out_dir <- file.path("benchmarks", "results", "marginsvm_speed_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workers <- max(1L, as.integer(Sys.getenv("SPATIAL_REFINE_BENCHMARK_WORKERS", "4")))
patterns <- c("jagged_stripes", "wavy_layers", "rings", "spiral", "branching",
              "lobes", "islands", "disconnected", "thin_layers", "intermixed",
              "layers3d")

general <- crossing(
  pattern = patterns,
  density_profile = c("uniform", "moderate", "strong", "extreme", "hotspot"),
  noise_type = c("random", "boundary", "patch", "region"),
  replicate = 1:2,
  classes = 5L
)
many_class <- crossing(
  pattern = c("jagged_stripes", "wavy_layers", "rings", "spiral", "layers3d"),
  density_profile = c("uniform", "extreme", "hotspot"),
  noise_type = c("random", "boundary", "patch", "region"),
  replicate = 1:2,
  classes = 19L
)
design <- bind_rows(general, many_class) |>
  mutate(scenario = row_number(), dimensions = if_else(pattern == "layers3d", 3L, 2L))

settings <- list(
  "Previous default" = list(overlap = 0.40, landmarks = 64L,
                            epochs = 10L, tv_iterations = 40L),
  "Faster default" = list(overlap = 0.25, landmarks = 48L,
                          epochs = 8L, tv_iterations = 24L)
)

run_one <- function(index) {
  d <- design[index, ]
  sim <- simulate_spatial_domains(
    8000L, d$pattern, k = d$classes, dimensions = d$dimensions, noise = 0.25,
    noise_type = d$noise_type, density_profile = d$density_profile,
    seed = 4100000L + index
  )
  initially_wrong <- sim$labels != sim$truth
  boundary <- sim$boundary_proximity <= quantile(sim$boundary_proximity, 0.2)
  result <- bind_rows(lapply(seq_along(settings), function(j) {
    control <- utils::modifyList(
      list(workers = 1L, seed = 4200000L + index), settings[[j]])
    timing <- system.time(pred <- refine_spatial_svm(
      sim$xy, sim$labels, sim$samples, control = control))
    correct <- pred == sim$truth
    data.frame(
      setting = names(settings)[j], accuracy = mean(correct),
      boundary_accuracy = mean(correct[boundary]),
      correction_recall = mean(correct[initially_wrong]),
      damage_rate = mean(!correct[!initially_wrong]),
      seconds = unname(timing[["elapsed"]])
    )
  }))
  bind_cols(d[rep(1L, length(settings)), ], result)
}

message("Running ", nrow(design), " independent validation scenarios")
rows <- if (.Platform$OS.type != "windows" && workers > 1L) {
  parallel::mclapply(seq_len(nrow(design)), run_one, mc.cores = workers,
                     mc.preschedule = TRUE)
} else {
  lapply(seq_len(nrow(design)), run_one)
}
metrics <- bind_rows(rows)
write.csv(metrics, file.path(out_dir, "speed_validation_metrics.csv"), row.names = FALSE)

paired <- metrics |>
  select(scenario, classes, pattern, density_profile, noise_type, setting,
         accuracy, boundary_accuracy, correction_recall, damage_rate, seconds) |>
  pivot_wider(names_from = setting,
              values_from = c(accuracy, boundary_accuracy, correction_recall,
                              damage_rate, seconds)) |>
  mutate(
    accuracy_delta = `accuracy_Faster default` - `accuracy_Previous default`,
    boundary_delta = `boundary_accuracy_Faster default` - `boundary_accuracy_Previous default`,
    correction_delta = `correction_recall_Faster default` - `correction_recall_Previous default`,
    damage_delta = `damage_rate_Faster default` - `damage_rate_Previous default`,
    speedup = `seconds_Previous default` / `seconds_Faster default`
  )
write.csv(paired, file.path(out_dir, "speed_validation_paired.csv"), row.names = FALSE)

set.seed(4300001L)
summarise_group <- function(data, label) {
  boot <- replicate(10000L, mean(sample(data$accuracy_delta, replace = TRUE)))
  data.frame(
    subset = label, scenarios = nrow(data),
    accuracy_delta = mean(data$accuracy_delta),
    accuracy_lower_95 = quantile(boot, 0.025),
    accuracy_upper_95 = quantile(boot, 0.975),
    boundary_delta = mean(data$boundary_delta),
    correction_delta = mean(data$correction_delta),
    damage_delta = mean(data$damage_delta),
    median_speedup = median(data$speedup),
    mean_speedup = mean(data$`seconds_Previous default`) /
      mean(data$`seconds_Faster default`)
  )
}
summary <- bind_rows(
  summarise_group(paired, "All"),
  summarise_group(filter(paired, classes == 5L), "5 classes"),
  summarise_group(filter(paired, classes == 19L), "19 classes")
)
write.csv(summary, file.path(out_dir, "speed_validation_summary.csv"), row.names = FALSE)

plot_data <- paired |>
  mutate(classes = paste(classes, "classes"))
plot <- ggplot(plot_data, aes(accuracy_delta, fill = classes)) +
  geom_histogram(position = "identity", alpha = 0.55, bins = 35) +
  geom_vline(xintercept = 0, color = "grey35", linewidth = 0.5) +
  facet_wrap(~classes, scales = "free_y") +
  labs(x = "Accuracy difference: faster minus previous default", y = "Scenarios", fill = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")
ggsave(file.path(out_dir, "validation_accuracy_differences.png"), plot,
       width = 8.2, height = 4.2, dpi = 240, bg = "white")

print(summary)
