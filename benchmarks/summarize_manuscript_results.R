#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

out_dir <- file.path("benchmarks", "results", "reviewer_response")
robust <- read.csv(file.path(out_dir, "robustness_metrics.csv"))
confirmatory <- robust |> filter(replicate > 3L)
write.csv(confirmatory, file.path(out_dir, "confirmatory_heldout_metrics.csv"), row.names = FALSE)

overall <- confirmatory |>
  group_by(method) |>
  summarise(
    accuracy = mean(accuracy), ari = mean(ari), boundary_010 = mean(boundary_010),
    smallest_recall = mean(smallest_recall), damage_rate = mean(damage_rate),
    seconds = mean(seconds), .groups = "drop"
  )
write.csv(overall, file.path(out_dir, "table_overall_confirmatory.csv"), row.names = FALSE)

by_error <- confirmatory |>
  filter(method %in% c("Initial", "Spatial graph", "SpaGCN refine", "GraphST refine",
                       "C++ kNN vote", "Potts-like ICM")) |>
  group_by(noise_type, method) |>
  summarise(
    accuracy = mean(accuracy), boundary_010 = mean(boundary_010),
    changed_precision = mean(changed_precision, na.rm = TRUE),
    correction_recall = mean(correction_recall, na.rm = TRUE),
    damage_rate = mean(damage_rate), .groups = "drop"
  )
write.csv(by_error, file.path(out_dir, "table_error_mechanisms.csv"), row.names = FALSE)

keys <- c("pattern", "noise", "noise_type", "replicate")
wide <- confirmatory |>
  select(all_of(keys), method, accuracy, boundary_010, damage_rate) |>
  pivot_wider(names_from = method, values_from = c(accuracy, boundary_010, damage_rate))
comparators <- c("Initial", "SpaGCN refine", "GraphST refine",
                 "C++ kNN vote", "Potts-like ICM")
set.seed(1500001L)
intervals <- list()
counter <- 0L
for (metric in c("accuracy", "boundary_010", "damage_rate")) {
  graph <- wide[[paste0(metric, "_Spatial graph")]]
  for (comparator in comparators) {
    difference <- graph - wide[[paste0(metric, "_", comparator)]]
    bootstrap <- replicate(10000L, mean(sample(difference, replace = TRUE), na.rm = TRUE))
    counter <- counter + 1L
    intervals[[counter]] <- data.frame(
      metric = metric, comparator = comparator,
      mean_difference = mean(difference, na.rm = TRUE),
      lower_95 = unname(quantile(bootstrap, 0.025, na.rm = TRUE)),
      upper_95 = unname(quantile(bootstrap, 0.975, na.rm = TRUE))
    )
  }
}
intervals <- bind_rows(intervals)
write.csv(intervals, file.path(out_dir, "table_paired_bootstrap_confirmatory.csv"), row.names = FALSE)

direct_intervals <- intervals |>
  filter(comparator != "Initial") |>
  mutate(
    favorable_difference = if_else(metric == "damage_rate", -mean_difference, mean_difference),
    favorable_lower = if_else(metric == "damage_rate", -upper_95, lower_95),
    favorable_upper = if_else(metric == "damage_rate", -lower_95, upper_95),
    metric = factor(
      metric,
      levels = c("accuracy", "boundary_010", "damage_rate"),
      labels = c("Overall accuracy", "Boundary accuracy", "Correct-label preservation")
    ),
    comparator = factor(
      comparator,
      levels = c("GraphST refine", "SpaGCN refine", "C++ kNN vote", "Potts-like ICM")
    )
  )
write.csv(direct_intervals, file.path(out_dir, "table_direct_paired_effects.csv"), row.names = FALSE)

p_direct <- ggplot(
  direct_intervals,
  aes(favorable_difference, comparator, xmin = favorable_lower, xmax = favorable_upper)
) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_errorbar(width = 0.18, linewidth = 0.55, orientation = "y") +
  geom_point(size = 2.2) +
  facet_wrap(~metric, scales = "free_x") +
  labs(
    x = "Paired difference (positive favors Spatial graph)", y = NULL,
    caption = "Points are mean paired differences; bars are 95% bootstrap intervals over matched held-out simulations."
  ) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey94"))
ggsave(file.path(out_dir, "direct_comparator_effects.png"), p_direct,
       width = 10, height = 4.2, dpi = 220)

support_summary <- confirmatory |>
  filter(method == "Spatial graph") |>
  summarise(
    expected_calibration_error = mean(ece, na.rm = TRUE),
    brier = mean(brier, na.rm = TRUE),
    changed_precision = mean(changed_precision, na.rm = TRUE),
    correction_recall = mean(correction_recall, na.rm = TRUE),
    damage_rate = mean(damage_rate)
  )
write.csv(support_summary, file.path(out_dir, "table_support_diagnostics.csv"), row.names = FALSE)

plot_summary <- confirmatory |>
  mutate(method = factor(method, c("Initial", "Spatial graph", "SpaGCN refine",
                                   "GraphST refine", "C++ kNN vote", "Potts-like ICM"))) |>
  group_by(noise_type, noise, method) |>
  summarise(accuracy = mean(accuracy), boundary_010 = mean(boundary_010), .groups = "drop")

p <- ggplot(plot_summary, aes(noise, accuracy, color = method)) +
  geom_line() + geom_point(size = 1.4) + facet_wrap(~noise_type) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Corrupted-label fraction", y = "Held-out accuracy", color = "Method") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "confirmatory_robustness.png"), p, width = 10, height = 5.5, dpi = 220)

p_boundary <- ggplot(plot_summary, aes(noise, boundary_010, color = method)) +
  geom_line() + geom_point(size = 1.4) + facet_wrap(~noise_type) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Corrupted-label fraction", y = "Held-out analytic-boundary accuracy", color = "Method") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "confirmatory_boundary.png"), p_boundary, width = 10, height = 5.5, dpi = 220)

message("Wrote held-out manuscript summaries")
