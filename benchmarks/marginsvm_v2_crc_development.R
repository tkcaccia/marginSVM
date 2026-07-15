#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SpatialGraphRefine)
  library(FNN)
  library(mclust)
  library(ggplot2)
})

object_path <- Sys.getenv("SPATIAL_REFINE_CRC_RDS",
  "/Users/stefano/Documents/wsitools/Data/VisiumHD/Colorectal/SpatialPolygons_2000genes_counts_only_with_tissue_annotations.rds")
change_threshold <- as.numeric(Sys.getenv("MARGINSVM_V2_CHANGE_THRESHOLD", "0.06"))
anisotropy <- as.numeric(Sys.getenv("MARGINSVM_V2_ANISOTROPY", "0.25"))
pairwise_specialists <- as.integer(Sys.getenv("MARGINSVM_V2_PAIRWISE", "1"))
threshold_tag <- gsub("\\.", "p", format(change_threshold, trim = TRUE))
phase <- Sys.getenv("MARGINSVM_V2_PHASE", "development")
seed_base <- as.integer(Sys.getenv("MARGINSVM_V2_SEED_BASE", "840000"))
out_dir <- file.path("benchmarks", "results",
  paste0("marginsvm_v2_crc_", phase, "_t", threshold_tag))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workers <- max(1L, as.integer(Sys.getenv("SPATIAL_REFINE_BENCHMARK_WORKERS", "4")))

object <- readRDS(object_path)
xy_all <- as.matrix(GetTissueCoordinates(object)[, c("x", "y")])
labels_all <- object@meta.data$wsi_annotation_name
keep <- !is.na(labels_all) & is.finite(xy_all[, 1L]) & is.finite(xy_all[, 2L])
xy <- xy_all[keep, , drop = FALSE]
truth <- droplevels(factor(labels_all[keep]))
code <- as.integer(truth)
n <- nrow(xy)
classes <- nlevels(truth)
samples <- factor(rep.int("CRC", n))

nn <- FNN::get.knn(xy, k = 8L)$nn.index
neighbor_code <- matrix(code[nn], nrow = n)
boundary <- rowSums(neighbor_code != code) > 0L
adjacent <- code
for (j in seq_len(ncol(nn))) {
  candidate <- neighbor_code[, j]
  use <- adjacent == code & candidate != code
  adjacent[use] <- candidate[use]
}
adjacent[adjacent == code] <- code[adjacent == code] %% classes + 1L
unit_xy <- sweep(xy, 2L, apply(xy, 2L, min), "-")
unit_xy <- sweep(unit_xy, 2L, apply(unit_xy, 2L, max), "/")

corrupt <- function(fraction, mechanism, seed) {
  set.seed(seed)
  count <- round(fraction * n)
  if (mechanism == "random") {
    index <- sample.int(n, count)
    replacement <- sample.int(classes - 1L, count, replace = TRUE)
    replacement <- replacement + (replacement >= code[index])
  } else if (mechanism == "boundary") {
    candidates <- which(boundary)
    index <- if (length(candidates) >= count) sample(candidates, count) else
      c(candidates, sample(setdiff(seq_len(n), candidates), count - length(candidates)))
    replacement <- adjacent[index]
  } else {
    centers <- sample.int(n, if (mechanism == "patch") 16L else 1L)
    distance <- vapply(centers, function(i) rowSums((unit_xy - unit_xy[i, ])^2), numeric(n))
    if (is.null(dim(distance))) distance <- matrix(distance, ncol = 1L)
    owner <- max.col(-distance)
    index <- order(distance[cbind(seq_len(n), owner)])[seq_len(count)]
    replacement <- adjacent[centers][owner[index]]
    same <- replacement == code[index]
    replacement[same] <- adjacent[index[same]]
  }
  initial <- truth
  initial[index] <- levels(truth)[replacement]
  list(labels = initial, corrupted = seq_len(n) %in% index)
}

design <- expand.grid(noise = c(0.10, 0.25),
  mechanism = c("random", "boundary", "patch", "region"), stringsAsFactors = FALSE)
scenario_filter <- as.integer(Sys.getenv("MARGINSVM_V2_SCENARIO", "0"))
scenario_ids <- if (scenario_filter %in% seq_len(nrow(design))) scenario_filter else seq_len(nrow(design))
metrics <- list()
predictions <- list()
position <- 0L
for (i in scenario_ids) {
  damaged <- corrupt(design$noise[i], design$mechanism[i], seed_base + i)
  for (method in c("marginSVM v1", "marginSVM v2")) {
    timing <- system.time(pred <- refine_spatial_svm(
      xy, damaged$labels, samples,
      control = list(experimental_v2 = as.integer(method == "marginSVM v2"),
                     change_threshold = change_threshold,
                     unresolved_threshold = min(0.04, change_threshold),
                     anisotropy = anisotropy,
                     pairwise_specialists = pairwise_specialists,
                     workers = workers, seed = seed_base + 10000L + i)
    ))
    correct <- pred == truth
    decision <- attr(pred, "decision")
    trust <- attr(pred, "trust")
    recalls <- vapply(levels(truth), function(label) mean(pred[truth == label] == label), numeric(1L))
    position <- position + 1L
    metrics[[position]] <- data.frame(
      scenario = i, noise = design$noise[i], mechanism = design$mechanism[i], method = method,
      accuracy = mean(correct), ari = adjustedRandIndex(pred, truth),
      macro_recall = mean(recalls), worst_recall = min(recalls),
      boundary_accuracy = mean(correct[boundary]), interior_accuracy = mean(correct[!boundary]),
      correction_recall = mean(correct[damaged$corrupted]),
      damage_rate = mean(!correct[!damaged$corrupted]), seconds = unname(timing["elapsed"]),
      mean_trust = if (is.null(trust)) 0 else mean(trust),
      unresolved = if (is.null(decision)) 0 else mean(decision == "unresolved"),
      stringsAsFactors = FALSE)
    predictions[[paste(i, method)]] <- pred
  }
  predictions[[paste(i, "Initial")]] <- damaged$labels
  message("completed CRC scenario ", i, "/", nrow(design))
}
metrics <- do.call(rbind, metrics)
write.csv(metrics, file.path(out_dir, "crc_development_metrics.csv"), row.names = FALSE)
summary <- aggregate(cbind(accuracy, ari, macro_recall, worst_recall, boundary_accuracy,
                           interior_accuracy, correction_recall, damage_rate, seconds,
                           mean_trust, unresolved) ~ method,
                     metrics, mean)
write.csv(summary, file.path(out_dir, "crc_development_summary.csv"), row.names = FALSE)

wide <- reshape(metrics[, c("scenario", "method", "accuracy")], idvar = "scenario",
                timevar = "method", direction = "wide")
best <- wide$scenario[which.max(wide[["accuracy.marginSVM v2"]] - wide[["accuracy.marginSVM v1"]])]
center_candidates <- which(boundary &
  predictions[[paste(best, "marginSVM v2")]] == truth &
  predictions[[paste(best, "marginSVM v1")]] != truth)
center <- if (length(center_candidates)) center_candidates[ceiling(length(center_candidates) / 2)] else which(boundary)[1L]
span <- c(diff(range(xy[, 1L])), diff(range(xy[, 2L]))) * 0.12
zoom <- abs(xy[, 1L] - xy[center, 1L]) <= span[1L] &
        abs(xy[, 2L] - xy[center, 2L]) <= span[2L]
if (sum(zoom) < 200L) zoom <- rep(TRUE, n)
display <- which(zoom)
if (length(display) > 35000L) display <- sample(display, 35000L)
panels <- c("Reference", "Corrupted input", "marginSVM v1", "marginSVM v2")
values <- list(truth, predictions[[paste(best, "Initial")]],
               predictions[[paste(best, "marginSVM v1")]],
               predictions[[paste(best, "marginSVM v2")]])
plot_data <- do.call(rbind, lapply(seq_along(panels), function(j) data.frame(
  x = xy[display, 1L], y = xy[display, 2L], label = values[[j]][display], panel = panels[j])))
plot_data$panel <- factor(plot_data$panel, panels)
p <- ggplot(plot_data, aes(x, y, color = label)) + geom_point(size = 0.28, alpha = 0.85) +
  facet_wrap(~panel, nrow = 1L) + coord_equal() + labs(x = NULL, y = NULL, color = "Annotation") +
  theme_void(base_size = 10) + theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = NA))
ggsave(file.path(out_dir, "crc_best_zoom.png"), p, width = 14, height = 5.2, dpi = 250, bg = "white")
print(summary)
