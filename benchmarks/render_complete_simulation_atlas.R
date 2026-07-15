#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(marginSVM)
  library(dplyr)
  library(ggplot2)
})

out_dir <- file.path("benchmarks", "results", "marginsvm_complete_simulation")
patterns <- c("jagged_stripes", "wavy_layers", "rings", "spiral", "branching",
              "lobes", "islands", "disconnected", "thin_layers", "intermixed",
              "layers3d", "gradient_A-B-B-C")

rows <- list()
position <- 0L
for (i in seq_along(patterns)) {
  pattern <- patterns[i]
  sim <- if (pattern == "gradient_A-B-B-C") {
    simulate_gradient_regions(6000, minority = 0.25, density_profile = "strong",
                              seed = 1600000L + i)
  } else {
    simulate_spatial_domains(
      6000, pattern, dimensions = if (pattern == "layers3d") 3L else 2L,
      noise = 0.25, noise_type = "random", density_profile = "strong",
      seed = 1600000L + i
    )
  }
  pred <- marginSVM:::.refine_spatial_svm_engine(
    sim$xy, sim$labels, sim$samples,
    control = list(workers = 1L, seed = 1700000L + i))
  display <- sample.int(nrow(sim$xy), 5000L)
  states <- list(Reference = sim$truth, `Corrupted input` = sim$labels,
                 marginSVM = pred)
  for (state in names(states)) {
    position <- position + 1L
    rows[[position]] <- data.frame(
      x = sim$xy[display, 1L], y = sim$xy[display, 2L],
      pattern = pattern, state = state, label = states[[state]][display]
    )
  }
}

atlas <- bind_rows(rows)
atlas$pattern <- factor(atlas$pattern, levels = patterns)
atlas$state <- factor(atlas$state, c("Reference", "Corrupted input", "marginSVM"))

for (part in seq_len(2L)) {
  selected <- patterns[((part - 1L) * 6L + 1L):(part * 6L)]
  plot_data <- filter(atlas, pattern %in% selected)
  plot_data$pattern <- droplevels(plot_data$pattern)
  plot <- ggplot(plot_data, aes(x, y, color = label)) +
    geom_point(size = 0.18, alpha = 0.82) + coord_equal() +
    facet_grid(pattern ~ state) + guides(color = "none") +
    labs(x = NULL, y = NULL) + theme_void(base_size = 9) +
    theme(strip.text.y = element_text(angle = 0),
          strip.text.x = element_text(face = "bold"),
          plot.background = element_rect(fill = "white", color = NA))
  ggsave(file.path(out_dir, paste0("complete_geometry_atlas_part", part, ".png")),
         plot, width = 9.5, height = 10.8, dpi = 240, bg = "white")
}
