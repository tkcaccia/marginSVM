#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

out_dir <- file.path("benchmarks", "results", "marginsvm_complete_simulation")
metrics <- read.csv(file.path(out_dir, "complete_simulation_metrics.csv"),
                    check.names = FALSE)
methods <- c("marginSVM", "Graph refinement",
             "SpaGCN refine", "GraphST refine", "C++ kNN vote",
             "Potts-like ICM")

# Replicates share the same simulation condition. Average them before inference
# so that seeds do not masquerade as independent benchmark tasks.
conditions <- metrics |>
  filter(method %in% methods) |>
  mutate(condition = if_else(
    family == "geometric",
    paste(family, pattern, density_profile, noise_type, sep = "::"),
    paste(family, dimensions, density_profile, minority, sep = "::")
  )) |>
  group_by(condition, family, pattern, density_profile, noise_type,
           dimensions, minority, method) |>
  summarise(accuracy = mean(accuracy), .groups = "drop")

wide <- conditions |>
  select(condition, family, method, accuracy) |>
  pivot_wider(names_from = method, values_from = accuracy)

set.seed(1900001L)
hierarchical_effects <- bind_rows(lapply(setdiff(methods, "marginSVM"), function(comparator) {
  difference <- wide[["marginSVM"]] - wide[[comparator]]
  bootstrap <- replicate(10000L, {
    sampled <- unlist(lapply(split(seq_len(nrow(wide)), wide$family), function(index) {
      sample(index, length(index), replace = TRUE)
    }), use.names = FALSE)
    mean(difference[sampled])
  })
  data.frame(
    comparator = comparator,
    conditions = length(difference),
    mean_difference = mean(difference),
    lower_95 = unname(quantile(bootstrap, 0.025)),
    upper_95 = unname(quantile(bootstrap, 0.975)),
    wins = mean(difference > 0),
    ties = mean(difference == 0),
    losses = mean(difference < 0)
  )
}))
write.csv(hierarchical_effects,
          file.path(out_dir, "condition_level_paired_effects.csv"), row.names = FALSE)

ranked <- conditions |>
  group_by(condition) |>
  mutate(rank = rank(-accuracy, ties.method = "average")) |>
  ungroup() |>
  group_by(method) |>
  summarise(mean_rank = mean(rank), median_rank = median(rank),
            first_fraction = mean(rank == 1), .groups = "drop") |>
  arrange(mean_rank)
write.csv(ranked, file.path(out_dir, "condition_level_ranks.csv"), row.names = FALSE)

ranked$method <- factor(ranked$method, levels = rev(ranked$method))
rank_plot <- ggplot(ranked, aes(mean_rank, method)) +
  geom_segment(aes(x = 1, xend = mean_rank, yend = method),
               color = "grey75", linewidth = 0.6) +
  geom_point(size = 2.5, color = "#C44E52") +
  geom_text(aes(label = sprintf("%.2f", mean_rank)), hjust = -0.45, size = 3.2) +
  scale_x_continuous(limits = c(1, 6), breaks = 1:6) +
  labs(x = "Mean rank across 236 replicate-averaged conditions", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
ggsave(file.path(out_dir, "condition_level_mean_rank.png"), rank_plot,
       width = 8.2, height = 4.6, dpi = 240, bg = "white")

accuracy_matrix <- wide[, methods]
friedman <- stats::friedman.test(as.matrix(accuracy_matrix))
pairwise <- stats::pairwise.wilcox.test(
  conditions$accuracy, conditions$method, paired = TRUE,
  p.adjust.method = "holm", exact = FALSE
)
test_output <- capture.output(
  list(friedman = friedman, pairwise_holm = pairwise)
)
while (length(test_output) > 0L && !nzchar(trimws(tail(test_output, 1L)))) {
  test_output <- head(test_output, -1L)
}
writeLines(sub("[[:space:]]+$", "", test_output),
           file.path(out_dir, "condition_level_tests.txt"))

by_density <- conditions |>
  group_by(density_profile, method) |>
  summarise(accuracy = mean(accuracy), .groups = "drop") |>
  group_by(density_profile) |>
  mutate(rank = rank(-accuracy, ties.method = "average")) |>
  ungroup()
write.csv(by_density, file.path(out_dir, "condition_level_density.csv"), row.names = FALSE)

print(hierarchical_effects)
print(ranked)
print(friedman)
