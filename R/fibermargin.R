#' Refine noisy categorical labels on a spatial domain
#'
#' `refine_spatial_labels()` applies FiberMargin to a coordinate-indexed label
#' field. The function is deterministic: no training phase is run and all route
#' constants are fixed internally. For multiclass data, each class receives a
#' two-sided local evidence score from short directional neighborhood sweeps; in
#' the binary case it applies the analogous dyadic-volume ballot rule.
#'
#' The method can process multiple specimens at once. Each specimen is refined
#' independently from the others; labels are never shared across specimens. A
#' specimen with no spatial extent is returned unchanged.
#'
#' @param xy Numeric matrix with two or three spatial coordinates per row.
#' @param labels Initial categorical assignment for every row of `xy`.
#' @param samples Optional specimen identifier. Different specimens are refined
#'   independently.
#' @param workers Optional CPU budget used across specimens and independent
#'   spatial charts. `NULL` uses up to four physical cores when possible.
#'
#' @return A factor with the levels and row names of `labels`. Attributes
#'   `candidate`, `margin_score`, `required`, `repair_margin`,
#'   `atlas_dispersion`, `isolation`, and `changed` contain pointwise
#'   diagnostics. `margin_score` is a local support contrast on the internal
#'   coordinate lattice, `repair_margin` is the difference between
#'   `margin_score` and the adaptive admission threshold, and `changed`
#'   marks updated sites. `isolation` is a deterministic local-gap
#'   protection factor in the multiclass route (one for the binary ballot). A
#'   candidate is accepted exactly when it differs from the observed label and
#'   `repair_margin` is nonnegative.
#' @export
#' @examples
#' sim <- simulate_gradient_regions(n = 2000, samples = 2)
#' refined <- refine_spatial_labels(sim$xy, sim$labels, sim$samples)
#' mean(refined == sim$truth)
refine_spatial_labels <- function(xy, labels, samples = NULL, workers = NULL) {
  .fiber_margin_engine(xy, labels, samples, workers = workers, control = list())
}

.fiber_margin_engine <- function(xy, labels, samples = NULL, workers = NULL,
                                 control = list()) {
  xy <- as.matrix(xy)
  if (!is.numeric(xy) || ncol(xy) < 2L || ncol(xy) > 3L ||
      nrow(xy) < 1L || any(!is.finite(xy))) {
    stop("`xy` must be a finite numeric matrix with two or three columns.",
         call. = FALSE)
  }
  storage.mode(xy) <- "double"
  if (length(labels) != nrow(xy) || anyNA(labels)) {
    stop("`labels` must contain one non-missing assignment per row of `xy`.",
         call. = FALSE)
  }
  labels <- as.factor(labels)
  if (is.null(samples)) samples <- rep.int(1L, nrow(xy))
  if (length(samples) != nrow(xy) || anyNA(samples)) {
    stop("`samples` must contain one non-missing identifier per row of `xy`.",
         call. = FALSE)
  }
  samples <- as.factor(samples)
  if (!is.list(control) || (length(control) && is.null(names(control)))) {
    stop("Internal FiberMargin control must be a named list.", call. = FALSE)
  }

  if (is.null(workers)) {
    cores <- parallel::detectCores(logical = FALSE)
    if (is.na(cores)) cores <- 1L
    workers <- max(1L, min(4L, as.integer(cores) - 1L))
  }
  workers <- as.integer(workers)
  if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("`workers` must be `NULL` or one positive integer.", call. = FALSE)
  }
  chart_cap <- if (ncol(xy) == 2L) 9L else 12L

  observed_samples <- droplevels(samples)
  if (nlevels(observed_samples) == 1L) {
    native_control <- control
    chart_workers <- min(workers, chart_cap)
    native_control$threads <- chart_workers
    result <- .Call(
      `_fibermargin_fiber_margin_cpp`,
      xy,
      as.integer(labels),
      rep.int(1L, nrow(xy)),
      native_control
    )
    refined <- factor(levels(labels)[result$labels], levels = levels(labels))
    names(refined) <- rownames(xy)
    attr(refined, "candidate") <- factor(
      levels(labels)[result$candidate], levels = levels(labels)
    )
    attr(refined, "margin_score") <- result$margin_score
    attr(refined, "required") <- result$required
    attr(refined, "repair_margin") <- result$margin_score - result$required
    attr(refined, "atlas_dispersion") <- result$atlas_dispersion
    attr(refined, "isolation") <- result$isolation
    attr(refined, "changed") <- result$changed
    attr(refined, "workers") <- chart_workers
    return(refined)
  }

  groups <- split(seq_len(nrow(xy)), observed_samples, drop = TRUE)
  process_workers <- if (.Platform$OS.type != "windows") {
    min(workers, length(groups))
  } else {
    1L
  }
  native_threads <- min(
    chart_cap, max(1L, workers %/% process_workers)
  )
  run_native <- function(rows) {
    native_control <- control
    native_control$threads <- native_threads
    .Call(
      `_fibermargin_fiber_margin_cpp`,
      xy[rows, , drop = FALSE],
      as.integer(labels[rows]),
      rep.int(1L, length(rows)),
      native_control
    )
  }
  if (process_workers > 1L) {
    results <- parallel::mclapply(
      groups, run_native, mc.cores = process_workers,
      mc.preschedule = TRUE, mc.set.seed = FALSE
    )
  } else {
    results <- lapply(groups, run_native)
  }

  refined_code <- as.integer(labels)
  candidate <- refined_code
  margin_score <- required <- atlas_dispersion <- rep.int(0, nrow(xy))
  isolation <- rep.int(1, nrow(xy))
  changed <- rep.int(FALSE, nrow(xy))
  for (index in seq_along(groups)) {
    rows <- groups[[index]]
    result <- results[[index]]
    refined_code[rows] <- result$labels
    candidate[rows] <- result$candidate
    margin_score[rows] <- result$margin_score
    required[rows] <- result$required
    atlas_dispersion[rows] <- result$atlas_dispersion
    isolation[rows] <- result$isolation
    changed[rows] <- result$changed
  }

  refined <- factor(levels(labels)[refined_code], levels = levels(labels))
  names(refined) <- rownames(xy)
  attr(refined, "candidate") <- factor(
    levels(labels)[candidate], levels = levels(labels)
  )
  attr(refined, "margin_score") <- margin_score
  attr(refined, "required") <- required
  attr(refined, "repair_margin") <- margin_score - required
  attr(refined, "atlas_dispersion") <- atlas_dispersion
  attr(refined, "isolation") <- isolation
  attr(refined, "changed") <- changed
  attr(refined, "workers") <- process_workers * native_threads
  refined
}
