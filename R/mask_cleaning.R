#' Clean a noisy categorical mask
#'
#' `clean_categorical_mask()` applies FiberMargin directly to a two-dimensional
#' pixel mask or three-dimensional voxel mask. Missing cells are ignored and
#' returned unchanged. The repair uses only voxel coordinates and label values
#' (no image intensities, model logits, or reference annotations).
#'
#' @param mask A categorical matrix or three-dimensional array. Numeric,
#'   character, logical, and factor masks are supported. Missing cells are
#'   treated as void locations.
#' @param samples Optional matrix or array with the same dimensions as `mask`.
#'   Different sample identifiers are cleaned independently.
#' @param workers Optional CPU budget used across samples and independent
#'   spatial charts. `NULL` uses up to four physical cores when possible.
#'
#' @return A cleaned mask with the dimensions, dimnames, storage type, and void
#'   locations of `mask`. Pointwise arrays are attached as attributes
#'   `candidate`, `margin_score`, `required`, `repair_margin`,
#'   `atlas_dispersion`, `isolation`, and `changed`.
#'   `margin_score` is a local support contrast on the mask lattice, while
#'   `repair_margin` is the evidence gap above the local acceptance threshold.
#' @export
#' @examples
#' mask <- matrix(c(rep("A", 50), rep("B", 50)), nrow = 10)
#' mask[5, 5] <- "B"
#' cleaned <- clean_categorical_mask(mask, workers = 1)
clean_categorical_mask <- function(mask, samples = NULL, workers = NULL) {
  dimensions <- dim(mask)
  if (is.null(dimensions) || !length(dimensions) %in% c(2L, 3L) ||
      any(dimensions < 1L)) {
    stop("`mask` must be a non-empty matrix or three-dimensional array.",
         call. = FALSE)
  }
  if (!is.atomic(mask) || is.complex(mask) || is.raw(mask)) {
    stop("`mask` must contain factor, character, logical, integer, or numeric labels.",
         call. = FALSE)
  }
  keep <- !is.na(mask)
  if (!any(keep)) {
    stop("`mask` must contain at least one non-missing label.", call. = FALSE)
  }

  sample_values <- NULL
  if (!is.null(samples)) {
    if (!identical(dim(samples), dimensions) || anyNA(samples[keep])) {
      stop("`samples` must match `mask` and be non-missing at every labelled cell.",
           call. = FALSE)
    }
    sample_values <- samples[keep]
  }

  coordinates <- arrayInd(which(keep), .dim = dimensions)
  colnames(coordinates) <- c("row", "column", "slice")[seq_len(ncol(coordinates))]
  refined <- refine_spatial_labels(
    xy = coordinates,
    labels = mask[keep],
    samples = sample_values,
    workers = workers
  )

  output <- .restore_mask_values(mask, keep, refined)
  attr(output, "candidate") <- .restore_mask_values(
    mask, keep, attr(refined, "candidate")
  )
  for (attribute in c(
    "margin_score", "required", "repair_margin", "atlas_dispersion", "isolation"
  )) {
    values <- array(NA_real_, dim = dimensions, dimnames = dimnames(mask))
    values[keep] <- attr(refined, attribute)
    attr(output, attribute) <- values
  }
  changed <- array(NA, dim = dimensions, dimnames = dimnames(mask))
  changed[keep] <- attr(refined, "changed")
  attr(output, "changed") <- changed
  attr(output, "workers") <- attr(refined, "workers")
  output
}

#' Generate a corrupted categorical mask
#'
#' Creates reproducible errors for evaluating post-hoc mask cleaning. The
#' mechanisms represent independent label errors, boundary-concentrated errors,
#' and one coherent overwritten patch.
#'
#' @param mask A reference categorical matrix or three-dimensional array.
#' @param mechanism One of `"impulse"`, `"boundary"`, or `"patch"`.
#' @param rate Fraction of non-missing cells whose labels are replaced.
#' @param seed One integer random seed.
#'
#' @return A mask with the same representation as `mask`. Exactly
#'   `round(rate * n)` non-missing cells are changed, subject to at least one and
#'   at most `n - 1` changes. One correct exemplar of every reference class is
#'   retained; a rate that makes this impossible is rejected.
#' @export
#' @examples
#' truth <- matrix(rep(c("A", "B"), each = 50), nrow = 10)
#' noisy <- corrupt_categorical_mask(truth, "boundary", rate = 0.15, seed = 4)
corrupt_categorical_mask <- function(
    mask, mechanism = c("impulse", "boundary", "patch"),
    rate = 0.15, seed = 1L) {
  mechanism <- match.arg(mechanism)
  dimensions <- dim(mask)
  if (is.null(dimensions) || !length(dimensions) %in% c(2L, 3L) ||
      any(dimensions < 1L)) {
    stop("`mask` must be a non-empty matrix or three-dimensional array.",
         call. = FALSE)
  }
  if (length(rate) != 1L || !is.numeric(rate) || !is.finite(rate) ||
      rate <= 0 || rate >= 1) {
    stop("`rate` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  seed <- as.integer(seed)
  if (length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be one finite integer.", call. = FALSE)
  }

  keep <- !is.na(mask)
  active <- which(keep)
  if (length(active) < 2L) {
    stop("`mask` must contain at least two non-missing cells.", call. = FALSE)
  }
  observed <- as.character(mask[active])
  classes <- unique(observed)
  if (length(classes) < 2L) {
    stop("`mask` must contain at least two observed classes.", call. = FALSE)
  }
  change_count <- min(
    length(active) - 1L,
    max(1L, as.integer(round(rate * length(active))))
  )
  set.seed(seed)
  active_position <- seq_along(active)
  boundary_mask <- as.logical(.categorical_mask_boundary(mask, keep))
  anchor_proximity <- as.numeric(!boundary_mask[active])

  if (mechanism == "impulse") {
    initial <- sample(active_position, change_count)
    selected_position <- .select_class_preserving_corruption(
      c(initial, setdiff(active_position, initial)), change_count, observed,
      boundary_proximity = anchor_proximity
    )
    selected <- active[selected_position]
    replacement <- vapply(as.character(mask[selected]), function(current) {
      sample(setdiff(classes, current), 1L)
    }, character(1L))
  } else if (mechanism == "boundary") {
    boundary_position <- which(boundary_mask[active])
    remainder <- setdiff(active_position, boundary_position)
    initial <- if (length(boundary_position) >= change_count) {
      sample(boundary_position, change_count)
    } else {
      c(boundary_position, sample(remainder, change_count - length(boundary_position)))
    }
    priority <- c(initial, setdiff(active_position, initial))
    selected_position <- .select_class_preserving_corruption(
      priority, change_count, observed, boundary_proximity = anchor_proximity
    )
    selected <- active[selected_position]
    replacement <- .adjacent_mask_replacements(mask, selected, keep, classes)
  } else {
    center <- sample(active, 1L)
    center_neighbors <- .adjacent_mask_replacements(
      mask, center, keep, classes
    )
    group <- interaction("sample1", observed, drop = TRUE, lex.order = TRUE)
    group_size <- table(group)
    target_candidates <- unique(c(center_neighbors, classes))
    feasible_target <- vapply(target_candidates, function(target) {
      candidate_group <- factor(group[observed != target], levels = levels(group))
      candidate_size <- table(candidate_group)
      capacity <- sum(pmin(as.integer(candidate_size), pmax(as.integer(group_size) - 1L, 0L)))
      capacity >= change_count
    }, logical(1L))
    if (!any(feasible_target)) {
      stop(
        "The requested corruption cannot retain one uncorrupted exemplar of every ",
        "true class within each specimen. Reduce the corruption rate.",
        call. = FALSE
      )
    }
    target <- target_candidates[[which(feasible_target)[1L]]]
    eligible_position <- which(observed != target)
    eligible <- active[eligible_position]
    center_coordinate <- arrayInd(center, .dim = dimensions)
    eligible_coordinates <- arrayInd(eligible, .dim = dimensions)
    squared_distance <- rowSums(
      (eligible_coordinates - matrix(
        center_coordinate, nrow(eligible_coordinates), length(dimensions),
        byrow = TRUE
      ))^2
    )
    priority <- eligible_position[order(squared_distance, runif(length(eligible)))]
    selected_position <- .select_class_preserving_corruption(
      priority, change_count, observed, boundary_proximity = anchor_proximity
    )
    selected <- active[selected_position]
    replacement <- rep(target, length(selected))
  }

  result <- mask
  source_values <- mask[active]
  source_keys <- as.character(source_values)
  mapped <- source_values[match(replacement, source_keys)]
  if (anyNA(mapped)) {
    stop("Internal corruption mapping failed.", call. = FALSE)
  }
  result[selected] <- mapped
  .validate_class_anchors(
    observed, result[active], context = "corrupt_categorical_mask()"
  )
  result
}

#' Evaluate categorical mask cleaning
#'
#' Evaluates a cleaned mask against a fixed reference and its imperfect input.
#' In addition to correction and damage metrics, the function reports overlap,
#' class-aware boundary agreement, and rare-class behavior.
#'
#' @param reference Reference categorical matrix or three-dimensional array.
#' @param initial Imperfect mask supplied to the cleaning method.
#' @param cleaned Output mask to evaluate.
#' @param elapsed Optional elapsed runtime in seconds.
#' @param method Optional method name included in the returned row.
#'
#' @details Missing cells in any input are excluded. Let `a0` be initial
#' accuracy, `r` correction recall, and `d` damage rate. The returned quantities
#' obey the exact decomposition
#' `accuracy = a0 * (1 - d) + (1 - a0) * r`.
#'
#' @return A one-row data frame of class-balanced, boundary, rare-class, repair,
#'   and damage measures.
#' @export
#' @examples
#' reference <- matrix(rep(c("A", "B"), each = 50), nrow = 10)
#' initial <- corrupt_categorical_mask(reference, seed = 8)
#' cleaned <- clean_categorical_mask(initial, workers = 1)
#' evaluate_mask_cleaning(reference, initial, cleaned)
evaluate_mask_cleaning <- function(reference, initial, cleaned,
                                   elapsed = NA_real_, method = NULL) {
  dimensions <- dim(reference)
  if (is.null(dimensions) || !length(dimensions) %in% c(2L, 3L) ||
      !identical(dim(initial), dimensions) ||
      !identical(dim(cleaned), dimensions)) {
    stop("`reference`, `initial`, and `cleaned` must have identical 2D or 3D dimensions.",
         call. = FALSE)
  }
  keep <- !is.na(reference) & !is.na(initial) & !is.na(cleaned)
  if (!any(keep)) {
    stop("The masks have no jointly observed cells.", call. = FALSE)
  }
  boundary <- .categorical_mask_boundary(reference, keep)[keep]
  reference_value <- as.character(reference[keep])
  initial_value <- as.character(initial[keep])
  cleaned_value <- as.character(cleaned[keep])
  class_count <- table(reference_value)
  rare_class <- names(class_count)[which.min(class_count)][1L]
  sparse <- reference_value == rare_class

  metrics <- evaluate_spatial_refinement(
    truth = reference_value,
    initial = initial_value,
    refined = cleaned_value,
    boundary = boundary,
    sparse = sparse,
    elapsed = elapsed,
    method = method
  )
  classes <- names(class_count)
  intersection <- vapply(classes, function(label) {
    sum(reference_value == label & cleaned_value == label)
  }, numeric(1L))
  union <- vapply(classes, function(label) {
    sum(reference_value == label | cleaned_value == label)
  }, numeric(1L))
  class_iou <- intersection / union
  class_dice <- 2 * intersection / vapply(classes, function(label) {
    sum(reference_value == label) + sum(cleaned_value == label)
  }, numeric(1L))
  class_boundary_iou <- vapply(classes, function(label) {
    reference_band <- .class_inner_boundary(reference, label, keep)
    cleaned_band <- .class_inner_boundary(cleaned, label, keep)
    band_union <- sum(reference_band | cleaned_band)
    if (band_union) {
      sum(reference_band & cleaned_band) / band_union
    } else {
      NA_real_
    }
  }, numeric(1L))
  repaired_fraction <- mean(
    initial_value != reference_value & cleaned_value == reference_value
  )
  damaged_fraction <- mean(
    initial_value == reference_value & cleaned_value != reference_value
  )

  extra <- data.frame(
    mean_iou = mean(class_iou),
    worst_iou = min(class_iou),
    frequency_weighted_iou = sum(class_iou * as.numeric(class_count)) /
      sum(class_count),
    macro_dice = mean(class_dice),
    mean_boundary_iou = if (all(is.na(class_boundary_iou))) {
      NA_real_
    } else {
      mean(class_boundary_iou, na.rm = TRUE)
    },
    rare_class = rare_class,
    rare_class_iou = unname(class_iou[rare_class]),
    boundary_fraction = mean(boundary),
    repaired_fraction = repaired_fraction,
    damaged_fraction = damaged_fraction,
    repair_identity_error = metrics$accuracy -
      (metrics$initial_accuracy * (1 - metrics$damage_rate) +
         (1 - metrics$initial_accuracy) * metrics$correction_recall),
    stringsAsFactors = FALSE
  )
  output <- cbind(metrics, extra)
  class(output) <- c("mask_cleaning_metrics", class(output))
  output
}

.restore_mask_values <- function(template, keep, values) {
  dimensions <- dim(template)
  dimension_names <- dimnames(template)
  if (is.factor(template)) {
    restored <- factor(
      rep.int(NA_character_, length(template)),
      levels = levels(template), ordered = is.ordered(template)
    )
    restored[keep] <- as.character(values)
    dim(restored) <- dimensions
    dimnames(restored) <- dimension_names
    return(restored)
  }

  restored <- template
  if (is.character(template)) {
    restored[keep] <- as.character(values)
  } else {
    source <- template[keep]
    mapped <- source[match(as.character(values), as.character(source))]
    if (anyNA(mapped)) {
      stop("A cleaned class could not be mapped to the input mask type.",
           call. = FALSE)
    }
    restored[keep] <- mapped
  }
  restored
}

.categorical_mask_boundary <- function(mask, keep = !is.na(mask)) {
  dimensions <- dim(mask)
  index <- array(seq_along(mask), dim = dimensions)
  values <- as.character(mask)
  boundary <- rep.int(FALSE, length(mask))
  for (axis in seq_along(dimensions)) {
    if (dimensions[axis] < 2L) next
    lower_subscript <- lapply(dimensions, seq_len)
    upper_subscript <- lower_subscript
    lower_subscript[[axis]] <- seq_len(dimensions[axis] - 1L)
    upper_subscript[[axis]] <- 2L:dimensions[axis]
    lower <- as.vector(do.call(`[`, c(list(index), lower_subscript, list(drop = FALSE))))
    upper <- as.vector(do.call(`[`, c(list(index), upper_subscript, list(drop = FALSE))))
    differs <- keep[lower] & keep[upper] & values[lower] != values[upper]
    boundary[c(lower[differs], upper[differs])] <- TRUE
  }
  array(boundary, dim = dimensions, dimnames = dimnames(mask))
}

.class_inner_boundary <- function(mask, label, keep = !is.na(mask)) {
  dimensions <- dim(mask)
  index <- array(seq_along(mask), dim = dimensions)
  inside <- array(
    keep & as.character(mask) == label,
    dim = dimensions, dimnames = dimnames(mask)
  )
  boundary <- rep.int(FALSE, length(mask))
  for (axis in seq_along(dimensions)) {
    if (dimensions[axis] < 2L) next
    lower_subscript <- lapply(dimensions, seq_len)
    upper_subscript <- lower_subscript
    lower_subscript[[axis]] <- seq_len(dimensions[axis] - 1L)
    upper_subscript[[axis]] <- 2L:dimensions[axis]
    lower <- as.vector(do.call(`[`, c(list(index), lower_subscript, list(drop = FALSE))))
    upper <- as.vector(do.call(`[`, c(list(index), upper_subscript, list(drop = FALSE))))
    differs <- inside[lower] != inside[upper]
    boundary[lower[differs & inside[lower]]] <- TRUE
    boundary[upper[differs & inside[upper]]] <- TRUE
  }
  array(boundary, dim = dimensions, dimnames = dimnames(mask))
}

.adjacent_mask_replacements <- function(mask, selected, keep, classes) {
  dimensions <- dim(mask)
  coordinates <- arrayInd(selected, .dim = dimensions)
  strides <- cumprod(c(1L, dimensions[-length(dimensions)]))
  values <- as.character(mask)
  vapply(seq_along(selected), function(index) {
    alternatives <- character()
    for (axis in seq_along(dimensions)) {
      if (coordinates[index, axis] > 1L) {
        neighbor <- selected[index] - strides[axis]
        if (keep[neighbor] && values[neighbor] != values[selected[index]]) {
          alternatives <- c(alternatives, values[neighbor])
        }
      }
      if (coordinates[index, axis] < dimensions[axis]) {
        neighbor <- selected[index] + strides[axis]
        if (keep[neighbor] && values[neighbor] != values[selected[index]]) {
          alternatives <- c(alternatives, values[neighbor])
        }
      }
    }
    alternatives <- unique(alternatives)
    if (!length(alternatives)) {
      alternatives <- setdiff(classes, values[selected[index]])
    }
    sample(alternatives, 1L)
  }, character(1L))
}
