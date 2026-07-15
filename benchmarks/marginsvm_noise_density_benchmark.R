#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(marginSVM)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

out_dir <- file.path("benchmarks", "results", "marginsvm_noise_density")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workers <- max(1L, as.integer(Sys.getenv("SPATIAL_REFINE_BENCHMARK_WORKERS", "4")))
patterns <- c("jagged_stripes", "wavy_layers", "rings", "spiral", "branching",
              "lobes", "islands", "disconnected", "thin_layers", "intermixed",
              "layers3d")
methods <- c("Initial", "marginSVM", "Graph refinement",
             "SpaGCN refine", "GraphST refine", "C++ kNN vote", "Potts-like ICM")
design <- tidyr::crossing(
  pattern = patterns,
  density_profile = c("uniform", "extreme"),
  noise_type = c("random", "boundary", "patch", "region"),
  noise = c(0.05, 0.15, 0.40),
  replicate = 1:3
) |>
  mutate(scenario = row_number(),
         dimensions = if_else(pattern == "layers3d", 3L, 2L),
         samples = if_else(replicate == 3L, 3L, 1L))

score <- function(pred, sim, method, seconds) {
  correct <- pred == sim$truth
  changed <- pred != sim$labels
  initially_wrong <- sim$labels != sim$truth
  counts <- table(droplevels(sim$region))
  sparse <- names(which.min(counts))
  data.frame(
    method = method,
    accuracy = mean(correct),
    sparse_region_accuracy = mean(correct[sim$region == sparse]),
    correction_recall = mean(correct[initially_wrong]),
    damage_rate = mean(!correct[!initially_wrong]),
    changed_precision = if (any(changed)) mean(correct[changed]) else NA_real_,
    seconds = seconds
  )
}

run_one <- function(index) {
  d <- design[index, ]
  sim <- simulate_spatial_domains(
    6000, d$pattern, dimensions = d$dimensions, samples = d$samples,
    noise = d$noise, noise_type = d$noise_type,
    density_profile = d$density_profile, seed = 2100000L + index
  )
  fns <- list(
    "marginSVM" = function() marginSVM:::.refine_spatial_svm_engine(
      sim$xy, sim$labels, sim$samples,
      control = list(workers = 1L, seed = 2200000L + index)),
    "Graph refinement" = function() marginSVM:::refine_spatial_clusters(
      sim$xy, sim$labels, sim$samples),
    "SpaGCN refine" = function() marginSVM:::.refine_published_labels(
      sim$xy, sim$labels, sim$samples, "spagcn", 6L),
    "GraphST refine" = function() marginSVM:::.refine_published_labels(
      sim$xy, sim$labels, sim$samples, "graphst", 50L),
    "C++ kNN vote" = function() marginSVM:::refine_spatial_clusters(
      sim$xy, sim$labels, sim$samples,
      control = list(weighted = FALSE, iterations = 1L, consensus = 0.5,
                     preserve = 0, margin = 0, current_support = 1)),
    "Potts-like ICM" = function() marginSVM:::refine_spatial_clusters(
      sim$xy, sim$labels, sim$samples,
      control = list(weighted = FALSE, iterations = 8L, consensus = 0.5,
                     preserve = 0.35, margin = 0, current_support = 1))
  )
  result <- list(score(sim$labels, sim, "Initial", 0))
  for (method in names(fns)) {
    timing <- system.time(pred <- fns[[method]]())
    result[[length(result) + 1L]] <- score(pred, sim, method, timing[["elapsed"]])
  }
  bind_cols(d[rep(1L, length(result)), ], bind_rows(result))
}

message("Running ", nrow(design), " noise-density scenarios")
results <- if (.Platform$OS.type != "windows" && workers > 1L) {
  parallel::mclapply(seq_len(nrow(design)), run_one, mc.cores = workers,
                     mc.preschedule = TRUE)
} else lapply(seq_len(nrow(design)), run_one)
new_metrics <- bind_rows(results)

complete <- read.csv(file.path(
  "benchmarks", "results", "marginsvm_complete_simulation",
  "complete_simulation_metrics.csv"), check.names = FALSE) |>
  filter(family == "geometric", density_profile %in% c("uniform", "extreme"),
         noise == 0.25) |>
  select(pattern, density_profile, noise_type, noise, replicate, dimensions,
         samples, method, accuracy, sparse_region_accuracy, correction_recall,
         damage_rate, changed_precision, seconds)
metrics <- bind_rows(
  new_metrics |> select(-scenario),
  complete
) |>
  mutate(method = factor(method, methods))
write.csv(metrics, file.path(out_dir, "noise_density_metrics.csv"), row.names = FALSE)

summary <- metrics |>
  filter(method != "Initial") |>
  group_by(density_profile, noise_type, noise, method) |>
  summarise(accuracy = mean(accuracy), sparse = mean(sparse_region_accuracy),
            damage = mean(damage_rate), .groups = "drop")
write.csv(summary, file.path(out_dir, "noise_density_summary.csv"), row.names = FALSE)

plot <- ggplot(summary, aes(noise, accuracy, color = method, group = method)) +
  geom_line(linewidth = 0.45) + geom_point(size = 1.0) +
  facet_grid(density_profile ~ noise_type) +
  scale_x_continuous(labels = scales::percent, breaks = c(0.05, 0.15, 0.25, 0.40)) +
  coord_cartesian(ylim = c(0.45, 1)) +
  labs(x = "Corrupted input labels", y = "Reference accuracy", color = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(out_dir, "accuracy_by_noise_density.png"), plot,
       width = 12.5, height = 6.4, dpi = 240, bg = "white")

print(summary |> filter(method == "marginSVM"))
