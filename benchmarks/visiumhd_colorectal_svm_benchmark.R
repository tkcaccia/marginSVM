#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(marginSVM)
  library(FNN)
  library(mclust)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

object_path <- Sys.getenv(
  "SPATIAL_REFINE_CRC_RDS",
  "/Users/stefano/Documents/wsitools/Data/VisiumHD/Colorectal/SpatialPolygons_2000genes_counts_only_with_tissue_annotations.rds"
)
phase <- match.arg(Sys.getenv("SPATIAL_REFINE_PHASE", "development"),
                   c("development", "heldout", "final"))
seed_base <- switch(phase, development = 240000L, heldout = 440000L,
                    final = 640000L)
out_dir <- file.path(
  "benchmarks", "results", paste0("visiumhd_colorectal_structured_svm_", phase)
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workers <- max(1L, as.integer(Sys.getenv("SPATIAL_REFINE_BENCHMARK_WORKERS", "4")))
quick <- identical(tolower(Sys.getenv("SPATIAL_REFINE_QUICK", "false")), "true")
figure_only <- identical(tolower(Sys.getenv("SPATIAL_REFINE_FIGURE_ONLY", "false")), "true")
target_tile_size <- max(1000L, as.integer(Sys.getenv("SPATIAL_REFINE_TARGET_TILE_SIZE", "10000")))
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank()))

message("Loading annotated VisiumHD colorectal section")
object <- readRDS(object_path)
xy_all <- as.matrix(GetTissueCoordinates(object)[, c("x", "y")])
labels_all <- object@meta.data$wsi_annotation_name
if (!identical(rownames(xy_all), rownames(object@meta.data))) {
  stop("Coordinate and metadata row orders do not match.")
}
keep <- !is.na(labels_all) & is.finite(xy_all[, 1L]) & is.finite(xy_all[, 2L])
xy <- xy_all[keep, , drop = FALSE]
truth <- droplevels(factor(labels_all[keep]))
truth_code <- as.integer(truth)
n <- nrow(xy)
n_classes <- nlevels(truth)
samples <- factor(rep.int("CRC_VisiumHD", n))

write.csv(
  data.frame(
    observations_total = nrow(xy_all), observations_annotated = n,
    observations_excluded = sum(!keep), classes = n_classes,
    duplicated_coordinates = sum(duplicated(xy))
  ),
  file.path(out_dir, "dataset_summary.csv"), row.names = FALSE
)
write.csv(
  data.frame(annotation = names(table(truth)), observations = as.integer(table(truth))),
  file.path(out_dir, "annotation_counts.csv"), row.names = FALSE
)

message("Building the boundary reference graph")
nn8 <- FNN::get.knn(xy, k = 8L)$nn.index
neighbor_codes <- matrix(truth_code[nn8], nrow = nrow(nn8))
boundary <- rowSums(neighbor_codes != truth_code) > 0L
adjacent_code <- truth_code
for (column in seq_len(ncol(nn8))) {
  candidate <- neighbor_codes[, column]
  replace <- adjacent_code == truth_code & candidate != truth_code
  adjacent_code[replace] <- candidate[replace]
}
adjacent_code[adjacent_code == truth_code] <-
  (truth_code[adjacent_code == truth_code] %% n_classes) + 1L

unit_xy <- sweep(xy, 2L, apply(xy, 2L, min), "-")
unit_xy <- sweep(unit_xy, 2L, pmax(apply(unit_xy, 2L, max), .Machine$double.eps), "/")

different_random_codes <- function(indices) {
  replacement <- sample.int(n_classes - 1L, length(indices), replace = TRUE)
  replacement + (replacement >= truth_code[indices])
}

corrupt_labels <- function(fraction, mechanism, seed) {
  set.seed(seed)
  n_noise <- min(n - 1L, as.integer(round(fraction * n)))
  if (mechanism == "random") {
    indices <- sample.int(n, n_noise)
    replacement <- different_random_codes(indices)
  } else if (mechanism == "boundary") {
    boundary_rows <- which(boundary)
    indices <- if (length(boundary_rows) >= n_noise) {
      sample(boundary_rows, n_noise)
    } else {
      c(boundary_rows, sample(setdiff(seq_len(n), boundary_rows), n_noise - length(boundary_rows)))
    }
    replacement <- adjacent_code[indices]
  } else if (mechanism == "patch") {
    centers <- sample.int(n, 16L)
    distance_matrix <- vapply(
      centers, function(center) rowSums((unit_xy - unit_xy[center, ])^2), numeric(n)
    )
    owner <- max.col(-distance_matrix, ties.method = "first")
    distance <- distance_matrix[cbind(seq_len(n), owner)]
    indices <- order(distance)[seq_len(n_noise)]
    replacement <- adjacent_code[centers][owner[indices]]
    same <- replacement == truth_code[indices]
    replacement[same] <- adjacent_code[indices[same]]
  } else if (mechanism == "region") {
    center <- sample.int(n, 1L)
    distance <- rowSums((unit_xy - unit_xy[center, ])^2)
    indices <- order(distance)[seq_len(n_noise)]
    replacement <- rep.int(adjacent_code[center], n_noise)
    same <- replacement == truth_code[indices]
    replacement[same] <- adjacent_code[indices[same]]
  } else {
    stop("Unknown corruption mechanism: ", mechanism)
  }
  output <- truth
  output[indices] <- levels(truth)[replacement]
  output
}

safe_mean <- function(x) if (length(x)) mean(x) else NA_real_
score_prediction <- function(pred, initial, method, elapsed) {
  correct <- pred == truth
  initially_wrong <- initial != truth
  initially_correct <- !initially_wrong
  changed <- pred != initial
  recalls <- vapply(levels(truth), function(annotation) {
    safe_mean(pred[truth == annotation] == annotation)
  }, numeric(1L))
  data.frame(
    method = method, accuracy = mean(correct),
    ari = mclust::adjustedRandIndex(pred, truth),
    macro_recall = mean(recalls), worst_recall = min(recalls),
    boundary_accuracy = safe_mean(correct[boundary]),
    interior_accuracy = safe_mean(correct[!boundary]),
    changed_precision = safe_mean(correct[changed]),
    correction_recall = safe_mean(correct[initially_wrong]),
    damage_rate = safe_mean(!correct[initially_correct]),
    changed_fraction = mean(changed), seconds = unname(elapsed),
    stringsAsFactors = FALSE
  )
}

# Frozen before this real-data evaluation from independent geometric simulations.
svm_settings <- list(
  target_tile_size = target_tile_size, workers = workers, seed = 230001L
)
svm_method <- function(initial) {
  marginSVM:::.refine_spatial_svm_engine(
    xy, initial, samples, control = svm_settings)
}
methods <- list(
  "SVM refinement" = svm_method,
  "Graph refinement" = function(initial) marginSVM:::refine_spatial_clusters(
    xy, initial, samples),
  "SpaGCN refine" = function(initial) marginSVM:::.refine_published_labels(
    xy, initial, samples, method = "spagcn", neighbors = 6L
  ),
  "GraphST refine" = function(initial) marginSVM:::.refine_published_labels(
    xy, initial, samples, method = "graphst", neighbors = 50L
  ),
  "C++ kNN vote" = function(initial) marginSVM:::refine_spatial_clusters(
    xy, initial, samples,
    control = list(weighted = FALSE, iterations = 1L, consensus = 0.5,
                   preserve = 0, margin = 0, current_support = 1)
  ),
  "Potts-like ICM" = function(initial) marginSVM:::refine_spatial_clusters(
    xy, initial, samples,
    control = list(weighted = FALSE, iterations = 8L, consensus = 0.5,
                   preserve = 0.35, margin = 0, current_support = 1)
  )
)

noise_levels <- if (quick) c(0.05, 0.25, 0.40) else c(0.05, 0.10, 0.15, 0.25, 0.40)
replicates <- if (quick) 1L else 1:3
test_design <- crossing(
  noise = noise_levels,
  mechanism = c("random", "boundary", "patch", "region"),
  replicate = replicates
)

message("Running frozen SVM and direct refiners on ", nrow(test_design), " corruptions")
run_test <- function(index, keep_predictions = FALSE) {
  scenario <- test_design[index, ]
  initial <- corrupt_labels(scenario$noise, scenario$mechanism, seed_base + index)
  rows <- list(score_prediction(initial, initial, "Initial", 0))
  predictions <- list(Initial = initial)
  for (method in names(methods)) {
    timing <- system.time(pred <- methods[[method]](initial))
    rows[[length(rows) + 1L]] <- score_prediction(pred, initial, method, timing["elapsed"])
    if (keep_predictions) predictions[[method]] <- pred
  }
  result <- bind_cols(scenario[rep(1L, length(rows)), ], bind_rows(rows))
  if (keep_predictions) attr(result, "predictions") <- predictions
  result
}

if (figure_only) {
  metrics <- read.csv(file.path(out_dir, "heldout_metrics.csv"))
} else {
  metrics <- bind_rows(lapply(seq_len(nrow(test_design)), function(index) {
    message("  scenario ", index, "/", nrow(test_design))
    run_test(index)
  }))
  write.csv(metrics, file.path(out_dir, "heldout_metrics.csv"), row.names = FALSE)
}

summary_metrics <- metrics |>
  group_by(method) |>
  summarise(
    accuracy = mean(accuracy), ari = mean(ari), macro_recall = mean(macro_recall),
    worst_recall = mean(worst_recall), boundary_accuracy = mean(boundary_accuracy),
    interior_accuracy = mean(interior_accuracy), correction_recall = mean(correction_recall),
    damage_rate = mean(damage_rate), changed_fraction = mean(changed_fraction),
    seconds = mean(seconds), .groups = "drop"
  ) |>
  arrange(desc(accuracy))
write.csv(summary_metrics, file.path(out_dir, "heldout_summary.csv"), row.names = FALSE)

paired <- metrics |>
  select(noise, mechanism, replicate, method, accuracy) |>
  pivot_wider(names_from = method, values_from = accuracy)
set.seed(250001L)
paired_intervals <- bind_rows(lapply(setdiff(names(methods), "SVM refinement"), function(comparator) {
  difference <- paired[["SVM refinement"]] - paired[[comparator]]
  bootstrap <- replicate(10000L, mean(sample(difference, replace = TRUE)))
  data.frame(
    comparator = comparator, mean_difference = mean(difference),
    lower_95 = unname(quantile(bootstrap, 0.025)),
    upper_95 = unname(quantile(bootstrap, 0.975))
  )
}))
write.csv(paired_intervals, file.path(out_dir, "svm_paired_effects.csv"), row.names = FALSE)

plot_metrics <- metrics |>
  group_by(noise, mechanism, method) |>
  summarise(accuracy = mean(accuracy), .groups = "drop")
p_accuracy <- ggplot(plot_metrics, aes(noise, accuracy, color = method)) +
  geom_line() + geom_point(size = 1.2) +
  facet_wrap(~mechanism, ncol = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Injected label-error fraction", y = "Reference annotation accuracy", color = "Method") +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "heldout_accuracy.png"), p_accuracy,
       width = 10, height = 7, dpi = 220, bg = "white")

p_tradeoff <- metrics |>
  filter(method != "Initial") |>
  group_by(method, mechanism, noise) |>
  summarise(
    correction_recall = mean(correction_recall), damage_rate = mean(damage_rate),
    .groups = "drop"
  ) |>
  ggplot(aes(damage_rate, correction_recall, color = method, shape = mechanism)) +
  geom_point(size = 2, alpha = 0.85) +
  labs(x = "Damage among initially correct labels", y = "Correction of injected errors",
       color = "Method", shape = "Error") +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "correction_damage_tradeoff.png"), p_tradeoff,
       width = 9, height = 6.5, dpi = 220, bg = "white")

example_index <- which(test_design$noise == 0.25 & test_design$mechanism == "random")[1L]
example_result <- run_test(example_index, keep_predictions = TRUE)
predictions <- attr(example_result, "predictions")
set.seed(260001L)
display <- sample.int(n, min(n, 70000L))
panels <- c("Reference" = "Reference", "Corrupted input" = "Initial",
            "marginSVM" = "SVM refinement", "GraphST refine" = "GraphST refine")
map_data <- bind_rows(lapply(names(panels), function(panel) {
  label <- if (panel == "Reference") truth else predictions[[panels[[panel]]]]
  data.frame(xy[display, , drop = FALSE], panel = panel, label = label[display])
}))
map_data$panel <- factor(map_data$panel, names(panels))
p_map <- ggplot(map_data, aes(x, y, color = label)) +
  geom_point(size = 0.12, alpha = 0.75) + facet_wrap(~panel, ncol = 2) +
  coord_equal() + labs(x = NULL, y = NULL, color = "Annotation") +
  theme_void(base_size = 9) +
  theme(
    legend.position = "bottom", legend.key.width = grid::unit(0.7, "lines"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
ggsave(file.path(out_dir, "random_25pct_spatial_example.png"), p_map,
       width = 10, height = 8, dpi = 240, bg = "white")

# Select zoom windows by frozen quantitative criteria rather than visual choice.
unit_display <- sweep(xy, 2L, apply(xy, 2L, min), "-")
unit_display <- sweep(unit_display, 2L, pmax(apply(unit_display, 2L, max),
                                            .Machine$double.eps), "/")
grid_x <- pmin(8L, floor(unit_display[, 1L] * 8L) + 1L)
grid_y <- pmin(6L, floor(unit_display[, 2L] * 6L) + 1L)
grid_id <- (grid_y - 1L) * 8L + grid_x
svm_correct <- predictions[["SVM refinement"]] == truth
graphst_correct <- predictions[["GraphST refine"]] == truth
initial_wrong <- predictions[["Initial"]] != truth
window_stats <- bind_rows(lapply(sort(unique(grid_id)), function(id) {
  rows <- grid_id == id
  composition <- prop.table(table(truth[rows]))
  data.frame(
    id = id, n = sum(rows),
    advantage = mean(svm_correct[rows]) - mean(graphst_correct[rows]),
    corrected = mean(initial_wrong[rows] & svm_correct[rows]),
    entropy = -sum(composition * log(pmax(composition, 1e-12)))
  )
})) |>
  filter(n >= 500L)
choose_distinct <- function(ordering, selected) {
  candidate <- ordering[!ordering %in% selected]
  candidate[1L]
}
selected <- integer()
selected <- c(selected, choose_distinct(window_stats$id[order(-window_stats$advantage)], selected))
selected <- c(selected, choose_distinct(window_stats$id[order(-window_stats$corrected)], selected))
selected <- c(selected, choose_distinct(window_stats$id[order(-window_stats$entropy)], selected))
zoom_names <- c("A: SVM advantage", "B: correction-rich", "C: multiclass border")
windows <- bind_rows(lapply(seq_along(selected), function(i) {
  rows <- grid_id == selected[i]
  padding_x <- 0.04 * diff(range(xy[, 1L]))
  padding_y <- 0.04 * diff(range(xy[, 2L]))
  data.frame(
    region = zoom_names[i],
    xmin = min(xy[rows, 1L]) - padding_x,
    xmax = max(xy[rows, 1L]) + padding_x,
    ymin = min(xy[rows, 2L]) - padding_y,
    ymax = max(xy[rows, 2L]) + padding_y
  )
}))
write.csv(windows, file.path(out_dir, "zoom_windows.csv"), row.names = FALSE)

overview_data <- data.frame(x = xy[display, 1L], y = xy[display, 2L],
                            label = truth[display])
p_overview <- ggplot(overview_data, aes(x, y, colour = label)) +
  geom_point(size = 0.10, alpha = 0.65) +
  geom_rect(
    data = windows,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, fill = NA, colour = "black", linewidth = 0.45
  ) +
  geom_text(
    data = windows,
    aes(x = xmin, y = ymax, label = substr(region, 1L, 1L)),
    inherit.aes = FALSE, hjust = -0.15, vjust = 1.15,
    colour = "black", fontface = "bold", size = 3.2
  ) +
  coord_equal() + theme_void(base_size = 9) + theme(legend.position = "none")

zoom_data <- bind_rows(lapply(seq_len(nrow(windows)), function(i) {
  box <- windows[i, ]
  rows <- xy[, 1L] >= box$xmin & xy[, 1L] <= box$xmax &
    xy[, 2L] >= box$ymin & xy[, 2L] <= box$ymax
  bind_rows(lapply(names(panels), function(panel) {
    label <- if (panel == "Reference") truth else predictions[[panels[[panel]]]]
    data.frame(x = xy[rows, 1L], y = xy[rows, 2L], label = label[rows],
               panel = panel, region = box$region)
  }))
}))
zoom_data$panel <- factor(zoom_data$panel, levels = names(panels))
zoom_data$region <- factor(zoom_data$region, levels = zoom_names)
p_zoom <- ggplot(zoom_data, aes(x, y, colour = label)) +
  geom_point(size = 0.34, alpha = 0.9) +
  coord_cartesian() +
  facet_wrap(vars(region, panel), ncol = 4, scales = "free") +
  theme_void(base_size = 9) +
  theme(
    strip.text = element_text(size = 8), legend.position = "bottom",
    aspect.ratio = 1,
    legend.key.width = grid::unit(0.6, "lines")
  ) +
  guides(colour = guide_legend(nrow = 3, override.aes = list(size = 2)))

p_zoom_combined <- p_overview / p_zoom + plot_layout(heights = c(0.65, 2.35))
ggsave(file.path(out_dir, "random_25pct_spatial_zoomed.png"), p_zoom_combined,
       width = 13.5, height = 11, dpi = 260, bg = "white")

message("VisiumHD colorectal SVM benchmark complete: ", out_dir)
