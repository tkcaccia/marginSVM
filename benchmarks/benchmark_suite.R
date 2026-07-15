#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialGraphRefine)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(FNN)
  library(mclust)
})

out_dir <- file.path("benchmarks", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_minimal(base_size = 12))

need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install benchmark dependency `", pkg, "` first.", call. = FALSE)
  }
}

majority <- function(x) {
  tab <- table(x)
  names(tab)[which.max(tab)]
}

spatial_knn_refine <- function(xy, labels, samples = NULL, k = 15L) {
  labels <- as.factor(labels)
  if (is.null(samples)) {
    samples <- factor(rep("sample1", length(labels)))
  } else {
    samples <- as.factor(samples)
  }
  pred <- labels
  for (s in levels(samples)) {
    idx <- which(samples == s)
    if (length(idx) <= 1L) next
    kk <- min(k + 1L, length(idx))
    nn <- FNN::get.knn(xy[idx, , drop = FALSE], k = kk)$nn.index
    pred[idx] <- apply(nn, 1L, function(ii) majority(labels[idx[ii]]))
  }
  factor(pred, levels = levels(labels))
}

e1071_svm_refine <- function(xy, labels, samples = NULL, gamma = 0.35, max_train = 6000L) {
  need_pkg("e1071")
  labels <- as.factor(labels)
  if (is.null(samples)) {
    samples <- factor(rep("sample1", length(labels)))
  } else {
    samples <- as.factor(samples)
  }
  pred <- labels
  for (s in levels(samples)) {
    idx <- which(samples == s)
    if (length(unique(labels[idx])) < 2L) next
    train_idx <- idx
    if (length(train_idx) > max_train) {
      train_idx <- sample(train_idx, max_train)
    }
    model <- e1071::svm(
      x = xy[train_idx, , drop = FALSE],
      y = labels[train_idx],
      kernel = "radial",
      gamma = gamma,
      scale = FALSE
    )
    pred[idx] <- as.character(predict(model, xy[idx, , drop = FALSE]))
  }
  factor(pred, levels = levels(labels))
}

spatial_continuity <- function(xy, labels, samples = NULL, k = 8L) {
  labels <- as.factor(labels)
  if (is.null(samples)) {
    samples <- factor(rep("sample1", length(labels)))
  } else {
    samples <- as.factor(samples)
  }
  vals <- numeric()
  for (s in levels(samples)) {
    idx <- which(samples == s)
    if (length(idx) <= 1L) next
    kk <- min(k + 1L, length(idx))
    nn <- FNN::get.knn(xy[idx, , drop = FALSE], k = kk)$nn.index
    local <- vapply(seq_len(nrow(nn)), function(i) {
      mean(labels[idx[nn[i, -1L]]] == labels[idx[i]])
    }, numeric(1))
    vals <- c(vals, local)
  }
  mean(vals)
}

simulate_domain_pattern <- function(n,
                                    dimensions = 2L,
                                    k = 5L,
                                    samples = 1L,
                                    noise = 0.1,
                                    pattern = c(
                                      "blobs", "unequal_blobs", "ellipses",
                                      "islands", "moons", "stripes",
                                      "jagged_stripes", "rings", "checker",
                                      "layers3d"
                                    ),
                                    seed = 1L) {
  pattern <- match.arg(pattern)
  if (pattern == "blobs") {
    return(simulate_spatial_clusters(n, dimensions, k, samples, noise, seed))
  }

  set.seed(seed)
  sample_id <- factor(rep(seq_len(samples), length.out = n))

  make_counts <- function(prob) {
    prob <- prob / sum(prob)
    counts <- as.integer(rmultinom(1L, n, prob))
    zero <- which(counts == 0L)
    if (length(zero) > 0L) {
      for (z in zero) {
        donor <- which.max(counts)
        counts[donor] <- counts[donor] - 1L
        counts[z] <- 1L
      }
    }
    counts
  }

  rotate2 <- function(x, angle) {
    rot <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), nrow = 2L)
    x %*% rot
  }

  if (pattern %in% c("unequal_blobs", "ellipses")) {
    angles <- seq(0, 2 * pi, length.out = k + 1L)[seq_len(k)]
    centers <- cbind(3.8 * cos(angles), 3.8 * sin(angles))
    if (pattern == "unequal_blobs") {
      probs <- rev(seq_len(k))^1.6
      counts <- make_counts(probs)
      sd_major <- seq(0.28, 0.95, length.out = k)
      sd_minor <- sd_major
    } else {
      counts <- make_counts(rep(1, k))
      sd_major <- seq(0.95, 1.45, length.out = k)
      sd_minor <- seq(0.12, 0.28, length.out = k)
    }
    xy2 <- matrix(NA_real_, nrow = n, ncol = 2L)
    truth_id <- integer(n)
    start <- 1L
    for (cl in seq_len(k)) {
      rows <- start:(start + counts[cl] - 1L)
      cloud <- cbind(rnorm(length(rows), sd = sd_major[cl]), rnorm(length(rows), sd = sd_minor[cl]))
      cloud <- rotate2(cloud, angles[cl] + pi / 5)
      xy2[rows, ] <- sweep(cloud, 2L, centers[cl, ], "+")
      truth_id[rows] <- cl
      start <- start + counts[cl]
    }
    ord <- sample.int(n)
    xy2 <- xy2[ord, , drop = FALSE]
    truth_id <- truth_id[ord]
    xy <- if (dimensions == 3L) cbind(xy2, rnorm(n, sd = 0.4)) else xy2
  } else if (pattern == "islands") {
    if (dimensions != 2L) stop("`islands` currently requires dimensions = 2.", call. = FALSE)
    probs <- c(0.58, rev(seq_len(k - 1L))^1.4)
    counts <- make_counts(probs)
    xy <- matrix(NA_real_, nrow = n, ncol = 2L)
    truth_id <- integer(n)
    rows <- seq_len(counts[1L])
    xy[rows, ] <- matrix(runif(length(rows) * 2L, min = -4, max = 4), ncol = 2L)
    truth_id[rows] <- 1L
    start <- counts[1L] + 1L
    island_centers <- cbind(runif(k - 1L, -3, 3), runif(k - 1L, -3, 3))
    for (cl in 2:k) {
      rows <- start:(start + counts[cl] - 1L)
      xy[rows, ] <- sweep(matrix(rnorm(length(rows) * 2L, sd = 0.22 + 0.04 * cl), ncol = 2L), 2L, island_centers[cl - 1L, ], "+")
      truth_id[rows] <- cl
      start <- start + counts[cl]
    }
    ord <- sample.int(n)
    xy <- xy[ord, , drop = FALSE]
    truth_id <- truth_id[ord]
  } else if (pattern == "moons") {
    if (dimensions != 2L || k != 2L) stop("`moons` requires dimensions = 2 and k = 2.", call. = FALSE)
    counts <- make_counts(c(1, 1))
    theta1 <- runif(counts[1L], 0, pi)
    theta2 <- runif(counts[2L], 0, pi)
    moon1 <- cbind(cos(theta1), sin(theta1)) + matrix(rnorm(counts[1L] * 2L, sd = 0.07), ncol = 2L)
    moon2 <- cbind(1 - cos(theta2), 0.45 - sin(theta2)) + matrix(rnorm(counts[2L] * 2L, sd = 0.07), ncol = 2L)
    xy <- rbind(moon1, moon2) * 3.2
    truth_id <- c(rep(1L, counts[1L]), rep(2L, counts[2L]))
    ord <- sample.int(n)
    xy <- xy[ord, , drop = FALSE]
    truth_id <- truth_id[ord]
  } else {
    xy <- matrix(runif(n * dimensions, min = -1, max = 1), ncol = dimensions)

    score <- switch(
      pattern,
      stripes = xy[, 1L] + 0.35 * sin(3 * pi * xy[, 2L]),
      jagged_stripes = xy[, 1L] + 0.28 * sin(5 * pi * xy[, 2L]) + 0.14 * sin(11 * pi * xy[, 2L]),
      rings = sqrt(rowSums(xy[, 1:2, drop = FALSE]^2)) + 0.08 * sin(8 * atan2(xy[, 2L], xy[, 1L])),
      checker = {
        gx <- floor((xy[, 1L] + 1) * k / 2)
        gy <- floor((xy[, 2L] + 1) * k / 2)
        (gx + gy) %% k
      },
      layers3d = {
        if (dimensions < 3L) stop("`layers3d` requires dimensions = 3.", call. = FALSE)
        xy[, 3L] + 0.25 * sin(2 * pi * xy[, 1L]) - 0.2 * cos(2 * pi * xy[, 2L])
      }
    )

    if (pattern == "checker") {
      truth_id <- as.integer(score) + 1L
    } else {
      cuts <- quantile(score, probs = seq(0, 1, length.out = k + 1L), names = FALSE)
      cuts <- unique(cuts)
      truth_id <- as.integer(cut(score, breaks = cuts, include.lowest = TRUE, labels = FALSE))
    }
    truth_id[is.na(truth_id)] <- 1L
  }

  k_actual <- max(truth_id)

  for (s in seq_len(samples)) {
    idx <- which(sample_id == s)
    xy[idx, ] <- sweep(xy[idx, , drop = FALSE], 2L, runif(dimensions, -0.15, 0.15), "+")
  }

  labels_id <- truth_id
  n_noise <- floor(noise * n)
  if (n_noise > 0L) {
    corrupt <- sample.int(n, n_noise)
    labels_id[corrupt] <- sample.int(k_actual, n_noise, replace = TRUE)
  }

  colnames(xy) <- paste0("dim", seq_len(dimensions))
  rownames(xy) <- paste0(pattern, seq_len(n))
  list(
    xy = xy,
    labels = factor(paste0("cluster", labels_id)),
    truth = factor(paste0("cluster", truth_id)),
    samples = sample_id
  )
}

metric_row <- function(method, pred, truth, xy, samples, elapsed, scenario) {
  scenario_name <- scenario$name
  scenario_k <- scenario$k
  scenario_samples <- scenario$samples
  scenario_noise <- scenario$noise
  scenario_tiles <- paste(scenario$tiles, collapse = "x")
  scenario_gamma <- scenario$gamma
  scenario_pattern <- scenario$pattern
  metric_accuracy <- mean(pred == truth)
  metric_ari <- mclust::adjustedRandIndex(pred, truth)
  metric_continuity <- spatial_continuity(xy, pred, samples)
  metric_seconds <- max(as.numeric(elapsed["elapsed"]), 1e-4)
  tibble(
    scenario = scenario_name,
    method = method,
    n = nrow(xy),
    dimensions = ncol(xy),
    clusters = scenario_k,
    samples = scenario_samples,
    pattern = scenario_pattern,
    noise = scenario_noise,
    tiles = scenario_tiles,
    gamma = scenario_gamma,
    seconds = metric_seconds,
    accuracy = metric_accuracy,
    ari = metric_ari,
    continuity = metric_continuity
  )
}

run_method <- function(label, expr) {
  gc()
  elapsed <- system.time(value <- force(expr))
  list(value = value, elapsed = elapsed)
}

plot_spatial_map <- function(sim, refined, scenario, suffix) {
  df <- as.data.frame(sim$xy)
  names(df) <- paste0("dim", seq_len(ncol(sim$xy)))
  df$truth <- sim$truth
  df$initial <- sim$labels
  df$refined <- refined
  if (ncol(sim$xy) == 3L) {
    df$z_bin <- cut(df$dim3, breaks = quantile(df$dim3, probs = seq(0, 1, length.out = 5)),
                    include.lowest = TRUE, labels = paste("z", 1:4))
  }
  long <- df |>
    select(dim1, dim2, truth, initial, refined, any_of("z_bin")) |>
    pivot_longer(c(truth, initial, refined), names_to = "state", values_to = "label")

  p <- ggplot(long, aes(dim1, dim2, color = label)) +
    geom_point(size = 0.25, alpha = 0.75) +
    coord_equal() +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(title = paste("Spatial maps:", scenario$name), x = "x", y = "y", color = "label")
  if ("z_bin" %in% names(long)) {
    p <- p + facet_grid(z_bin ~ state)
  } else {
    p <- p + facet_wrap(~state, nrow = 1L)
  }
  ggsave(file.path(out_dir, paste0("spatial_map_", suffix, ".png")), p, width = 12, height = if (ncol(sim$xy) == 3L) 10 else 4, dpi = 180)
}

scenarios <- tibble::tribble(
  ~name,              ~pattern,   ~n,     ~dimensions, ~k, ~samples, ~noise, ~tiles,        ~gamma, ~include_e1071,
  "2d_low_noise",     "blobs",     8000L,  2L,          5L, 1L,       0.05,   list(c(5, 5)), 0.35, TRUE,
  "2d_high_noise",    "blobs",     8000L,  2L,          5L, 1L,       0.25,   list(c(5, 5)), 0.35, TRUE,
  "2d_multi_sample",  "blobs",    12000L,  2L,          6L, 3L,       0.15,   list(c(6, 6)), 0.35, TRUE,
  "2d_many_domains",  "blobs",    14000L,  2L,         10L, 1L,       0.18,   list(c(8, 8)), 0.40, TRUE,
  "2d_few_domains",   "blobs",     9000L,  2L,          3L, 1L,       0.18,   list(c(4, 4)), 0.30, TRUE,
  "2d_imbalanced",    "unequal_blobs", 12000L, 2L,      7L, 1L,       0.20,   list(c(7, 7)), 0.40, TRUE,
  "2d_ellipses",      "ellipses", 12000L,  2L,          6L, 1L,       0.20,   list(c(6, 6)), 0.45, TRUE,
  "2d_rare_islands",  "islands",  12000L,  2L,          6L, 1L,       0.20,   list(c(7, 7)), 0.55, TRUE,
  "2d_moons",         "moons",    10000L,  2L,          2L, 1L,       0.18,   list(c(6, 6)), 0.60, TRUE,
  "3d_low_noise",     "blobs",     9000L,  3L,          5L, 1L,       0.08,   list(c(4, 4, 4)), 0.35, TRUE,
  "3d_high_noise",    "blobs",     9000L,  3L,          5L, 2L,       0.25,   list(c(4, 4, 4)), 0.35, TRUE,
  "2d_stripes",       "stripes",  10000L,  2L,          5L, 1L,       0.18,   list(c(6, 6)), 0.50, TRUE,
  "2d_jagged_stripes", "jagged_stripes", 10000L, 2L,    5L, 1L,       0.22,   list(c(8, 8)), 0.65, TRUE,
  "2d_rings",         "rings",    10000L,  2L,          5L, 1L,       0.18,   list(c(6, 6)), 0.50, TRUE,
  "2d_checker",       "checker",  10000L,  2L,          5L, 1L,       0.18,   list(c(8, 8)), 0.75, TRUE,
  "3d_layers",        "layers3d", 12000L,  3L,          5L, 2L,       0.18,   list(c(4, 4, 4)), 0.50, TRUE,
  "2d_60k_scale",     "blobs",    60000L,  2L,          6L, 2L,       0.12,   list(c(8, 8)), 0.35, FALSE,
  "3d_60k_scale",     "blobs",    60000L,  3L,          6L, 2L,       0.12,   list(c(5, 5, 5)), 0.35, FALSE
)

noise_sweep <- tidyr::expand_grid(
  pattern = c("blobs", "unequal_blobs", "ellipses", "islands", "jagged_stripes"),
  noise = c(0.05, 0.15, 0.30, 0.45)
) |>
  mutate(
    name = paste0("noise_", pattern, "_", sprintf("%02d", as.integer(noise * 100))),
    n = 9000L,
    dimensions = 2L,
    k = dplyr::case_when(
      pattern == "blobs" ~ 6L,
      pattern == "islands" ~ 6L,
      TRUE ~ 5L
    ),
    samples = 1L,
    tiles = dplyr::case_when(
      pattern == "jagged_stripes" ~ list(c(8L, 8L)),
      pattern == "islands" ~ list(c(7L, 7L)),
      TRUE ~ list(c(6L, 6L))
    ),
    gamma = dplyr::case_when(
      pattern == "jagged_stripes" ~ 0.65,
      pattern == "islands" ~ 0.55,
      TRUE ~ 0.45
    ),
    include_e1071 = FALSE
  ) |>
  select(name, pattern, n, dimensions, k, samples, noise, tiles, gamma, include_e1071)

size_sweep <- tidyr::expand_grid(
  k = c(3L, 5L, 8L, 12L),
  pattern = c("blobs", "unequal_blobs")
) |>
  mutate(
    name = paste0("size_", pattern, "_k", k),
    n = 10000L,
    dimensions = 2L,
    samples = 1L,
    noise = 0.20,
    tiles = lapply(k, function(kk) rep(max(4L, min(9L, kk)), 2L)),
    gamma = 0.40,
    include_e1071 = FALSE
  ) |>
  select(name, pattern, n, dimensions, k, samples, noise, tiles, gamma, include_e1071)

tile_sweep <- tidyr::expand_grid(
  tile_n = c(2L, 4L, 6L, 8L),
  gamma = c(0.1, 0.35, 1.0)
) |>
  mutate(
    name = paste0("tile", tile_n, "_gamma", gamma),
    n = 12000L,
    dimensions = 2L,
    k = 6L,
    samples = 2L,
    pattern = "blobs",
    noise = 0.18,
    tiles = lapply(tile_n, function(x) c(x, x)),
    include_e1071 = FALSE
  ) |>
  select(name, pattern, n, dimensions, k, samples, noise, tiles, gamma, include_e1071)

all_scenarios <- bind_rows(scenarios, noise_sweep, size_sweep, tile_sweep)

all_metrics <- list()

set.seed(20260630)
for (i in seq_len(nrow(all_scenarios))) {
  scenario <- as.list(all_scenarios[i, ])
  scenario$tiles <- unlist(scenario$tiles, use.names = FALSE)
  message("\nRunning ", scenario$name, " (n=", scenario$n, ", d=", scenario$dimensions, ")")
  sim <- simulate_domain_pattern(
    n = scenario$n,
    dimensions = scenario$dimensions,
    k = scenario$k,
    samples = scenario$samples,
    noise = scenario$noise,
    pattern = scenario$pattern,
    seed = 1000 + i
  )

  raw <- run_method("initial_labels", sim$labels)
  cpp <- run_method("SpatialGraphRefine_RFF", refine_spatial_svm(
    sim$xy,
    sim$labels,
    samples = sim$samples,
    tiles = scenario$tiles,
    gamma = scenario$gamma,
    n_features = 96,
    epochs = 6,
    seed = 42
  ))
  knn <- run_method("spatial_knn_majority", spatial_knn_refine(sim$xy, sim$labels, sim$samples, k = 15L))

  rows <- list(
    metric_row("initial_labels", raw$value, sim$truth, sim$xy, sim$samples, raw$elapsed, scenario),
    metric_row("SpatialGraphRefine_RFF", cpp$value, sim$truth, sim$xy, sim$samples, cpp$elapsed, scenario),
    metric_row("spatial_knn_majority", knn$value, sim$truth, sim$xy, sim$samples, knn$elapsed, scenario)
  )

  if (isTRUE(scenario$include_e1071)) {
    radial <- run_method("e1071_radial_svm_subsample", e1071_svm_refine(
      sim$xy, sim$labels, sim$samples, gamma = scenario$gamma, max_train = 6000L
    ))
    rows[[length(rows) + 1L]] <- metric_row(
      "e1071_radial_svm_subsample", radial$value, sim$truth, sim$xy, sim$samples, radial$elapsed, scenario
    )
  }

  all_metrics[[length(all_metrics) + 1L]] <- bind_rows(rows)

  if (scenario$name %in% c(
    "2d_high_noise", "2d_many_domains", "2d_imbalanced", "2d_ellipses",
    "2d_rare_islands", "2d_moons", "3d_high_noise", "2d_stripes",
    "2d_jagged_stripes", "2d_rings", "2d_checker", "3d_layers",
    "2d_60k_scale"
  )) {
    plot_spatial_map(sim, cpp$value, scenario, scenario$name)
  }
}

metrics <- bind_rows(all_metrics)
write.csv(metrics, file.path(out_dir, "benchmark_metrics.csv"), row.names = FALSE)

main_metrics <- metrics |>
  filter(!grepl("^tile", scenario))

p_accuracy <- ggplot(main_metrics, aes(method, accuracy, fill = method)) +
  geom_col(width = 0.75) +
  facet_wrap(~scenario, scales = "free_y") +
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none") +
  labs(title = "Cluster refinement accuracy across simulated spatial benchmarks", x = NULL, y = "accuracy")
ggsave(file.path(out_dir, "accuracy_by_scenario.png"), p_accuracy, width = 12, height = 7, dpi = 180)

p_ari <- ggplot(main_metrics, aes(method, ari, fill = method)) +
  geom_col(width = 0.75) +
  facet_wrap(~scenario, scales = "free_y") +
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none") +
  labs(title = "Adjusted Rand Index across simulated spatial benchmarks", x = NULL, y = "ARI")
ggsave(file.path(out_dir, "ari_by_scenario.png"), p_ari, width = 12, height = 7, dpi = 180)

p_runtime <- ggplot(main_metrics, aes(n, seconds, color = method)) +
  geom_point(size = 2) +
  geom_line(aes(group = interaction(method, dimensions))) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~dimensions, labeller = label_both) +
  labs(title = "Runtime scaling", x = "observations", y = "seconds (log scale)", color = "method")
ggsave(file.path(out_dir, "runtime_scaling.png"), p_runtime, width = 9, height = 5, dpi = 180)

p_continuity <- ggplot(main_metrics, aes(accuracy, continuity, color = method, size = seconds)) +
  geom_point(alpha = 0.85) +
  facet_wrap(~scenario) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Accuracy versus spatial continuity", x = "accuracy", y = "same-label neighbor fraction")
ggsave(file.path(out_dir, "accuracy_vs_continuity.png"), p_continuity, width = 11, height = 7, dpi = 180)

structure_metrics <- main_metrics |>
  filter(!grepl("^(noise_|size_)", scenario), !grepl("_60k_scale$", scenario))

p_structure <- ggplot(structure_metrics, aes(pattern, accuracy, fill = method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  facet_wrap(~scenario, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "Refinement performance across spatial structures", x = "structure", y = "accuracy", fill = "method")
ggsave(file.path(out_dir, "accuracy_by_structure.png"), p_structure, width = 13, height = 8, dpi = 180)

noise_metrics <- main_metrics |>
  filter(grepl("^noise_", scenario))

p_noise <- ggplot(noise_metrics, aes(noise, accuracy, color = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~pattern) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Robustness to increasing label noise", x = "label noise fraction", y = "accuracy", color = "method")
ggsave(file.path(out_dir, "noise_robustness.png"), p_noise, width = 11, height = 7, dpi = 180)

size_metrics <- main_metrics |>
  filter(grepl("^size_", scenario))

p_size <- ggplot(size_metrics, aes(clusters, accuracy, color = method, linetype = pattern)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Sensitivity to number of domains and domain-size imbalance", x = "number of domains", y = "accuracy", color = "method", linetype = "structure")
ggsave(file.path(out_dir, "cluster_count_size_sensitivity.png"), p_size, width = 10, height = 6, dpi = 180)

tile_metrics <- metrics |>
  filter(grepl("^tile", scenario), method == "SpatialGraphRefine_RFF") |>
  mutate(tile_n = as.integer(sub("x.*", "", tiles)))

p_tile <- ggplot(tile_metrics, aes(factor(tile_n), factor(gamma), fill = accuracy)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", accuracy)), size = 3) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(title = "Tile/gamma sensitivity for the RFF-SVM backend", x = "tiles per dimension", y = "gamma", fill = "accuracy")
ggsave(file.path(out_dir, "tile_gamma_heatmap.png"), p_tile, width = 7, height = 5, dpi = 180)

summary_table <- metrics |>
  arrange(scenario, desc(accuracy), seconds) |>
  group_by(scenario) |>
  slice_head(n = 4L) |>
  ungroup()
write.csv(summary_table, file.path(out_dir, "benchmark_summary_top_methods.csv"), row.names = FALSE)

message("\nWrote benchmark tables and figures to: ", normalizePath(out_dir))
