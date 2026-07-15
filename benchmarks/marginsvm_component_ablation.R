#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(marginSVM)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

out_dir <- file.path("benchmarks", "results", "marginsvm_component_ablation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workers <- max(1L, as.integer(Sys.getenv("SPATIAL_REFINE_BENCHMARK_WORKERS", "4")))
design <- crossing(
  pattern = c("jagged_stripes", "rings", "islands", "intermixed", "layers3d"),
  density_profile = c("uniform", "extreme"),
  noise_type = c("random", "boundary"),
  replicate = 1:3
) |>
  mutate(scenario = row_number(),
         dimensions = if_else(pattern == "layers3d", 3L, 2L))

variants <- list(
  "Full marginSVM" = list(),
  "No cross-fitting" = list(cross_fitting = 0),
  "Regular tiles only" = list(adaptive_tiles = 0),
  "No tile overlap" = list(overlap = 0),
  "No graph fusion" = list(graph_mix = 0),
  "No TV decoding" = list(tv_strength = 0),
  "No robust ramp" = list(ramp = 1e6)
)

run_one <- function(index) {
  d <- design[index, ]
  sim <- simulate_spatial_domains(
    6000, d$pattern, dimensions = d$dimensions, noise = 0.25,
    noise_type = d$noise_type, density_profile = d$density_profile,
    seed = 2300000L + index
  )
  rows <- bind_rows(lapply(seq_along(variants), function(v) {
    control <- utils::modifyList(
      list(workers = 1L, seed = 2400000L + index), variants[[v]])
    timing <- system.time(pred <- marginSVM:::.refine_spatial_svm_engine(
      sim$xy, sim$labels, sim$samples, control = control))
    correct <- pred == sim$truth
    initially_wrong <- sim$labels != sim$truth
    data.frame(
      variant = names(variants)[v],
      accuracy = mean(correct),
      boundary_accuracy = mean(correct[sim$boundary_proximity <=
        quantile(sim$boundary_proximity, 0.2)]),
      correction_recall = mean(correct[initially_wrong]),
      damage_rate = mean(!correct[!initially_wrong]),
      seconds = timing[["elapsed"]]
    )
  }))
  bind_cols(d[rep(1L, length(variants)), ], rows)
}

message("Running ", nrow(design), " ablation scenarios")
results <- if (.Platform$OS.type != "windows" && workers > 1L) {
  parallel::mclapply(seq_len(nrow(design)), run_one, mc.cores = workers,
                     mc.preschedule = TRUE)
} else lapply(seq_len(nrow(design)), run_one)
metrics <- bind_rows(results)
write.csv(metrics, file.path(out_dir, "component_ablation_metrics.csv"), row.names = FALSE)

summary <- metrics |>
  group_by(variant) |>
  summarise(accuracy = mean(accuracy), boundary_accuracy = mean(boundary_accuracy),
            correction_recall = mean(correction_recall), damage_rate = mean(damage_rate),
            seconds = mean(seconds), .groups = "drop") |>
  arrange(desc(accuracy))
write.csv(summary, file.path(out_dir, "component_ablation_summary.csv"), row.names = FALSE)

paired <- metrics |>
  select(scenario, variant, accuracy) |>
  pivot_wider(names_from = variant, values_from = accuracy)
set.seed(2401001L)
effects <- bind_rows(lapply(setdiff(names(variants), "Full marginSVM"), function(variant) {
  difference <- paired[["Full marginSVM"]] - paired[[variant]]
  bootstrap <- replicate(10000L, mean(sample(difference, replace = TRUE)))
  data.frame(variant = variant, difference = mean(difference),
             lower_95 = quantile(bootstrap, 0.025),
             upper_95 = quantile(bootstrap, 0.975))
}))
write.csv(effects, file.path(out_dir, "component_ablation_effects.csv"), row.names = FALSE)

plot <- summary |>
  mutate(variant = reorder(variant, accuracy)) |>
  ggplot(aes(accuracy, variant)) +
  geom_segment(aes(x = min(accuracy) - 0.005, xend = accuracy, yend = variant),
               color = "grey75", linewidth = 0.6) +
  geom_point(size = 2.4, color = "#C44E52") +
  labs(x = "Reference accuracy", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
ggsave(file.path(out_dir, "component_ablation_accuracy.png"), plot,
       width = 8.2, height = 4.8, dpi = 240, bg = "white")

print(summary)
