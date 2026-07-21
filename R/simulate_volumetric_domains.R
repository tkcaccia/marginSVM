#' Simulate complex three-dimensional tissue domains
#'
#' Generates volumetric spatial-label benchmarks with curved interfaces,
#' disconnected components, thin structures, class imbalance, irregular
#' z-plane acquisition, and spatially structured label corruption.
#'
#' @param n Number of observations.
#' @param shape Three-dimensional domain geometry.
#' @param acquisition Sampling design: `"uniform"`, `"class_imbalanced"`, or
#'   `"irregular_z"`.
#' @param noise_type Corruption mechanism: `"random"`, `"boundary"`,
#'   `"patch"`, or `"region"`.
#' @param noise Fraction of labels to corrupt within each sample. Each true
#'   class present in a specimen retains one correct exemplar; a rate that makes
#'   this impossible is rejected.
#' @param k Number of domains.
#' @param samples Number of independent tissue volumes.
#' @param seed Random seed.
#'
#' @return A `spatial_refinement_benchmark` with three coordinate columns.
#' @export
#' @examples
#' sim <- simulate_volumetric_domains(
#'   n = 1000, shape = "folded_layers", noise_type = "boundary", seed = 8
#' )
#' dim(sim$xy)
simulate_volumetric_domains <- function(
    n = 50000L,
    shape = c(
      "concentric_shells", "warped_ellipsoids", "toroidal_compartments",
      "folded_layers", "thin_folded_sheets", "branching_tubes",
      "disconnected_volumes", "helical_channels"
    ),
    acquisition = c("uniform", "class_imbalanced", "irregular_z"),
    noise_type = c("random", "boundary", "patch", "region"),
    noise = 0.25,
    k = 5L,
    samples = 1L,
    seed = 1L) {
  shape_names <- c(
    concentric_shells = "Concentric shells",
    warped_ellipsoids = "Warped ellipsoids",
    toroidal_compartments = "Toroidal compartments",
    folded_layers = "Folded layers",
    thin_folded_sheets = "Thin folded sheets",
    branching_tubes = "Branching tubes",
    disconnected_volumes = "Disconnected volumes",
    helical_channels = "Helical channels"
  )
  requested_shape <- shape[1L]
  shape_key <- if (requested_shape %in% unname(shape_names)) {
    names(shape_names)[match(requested_shape, shape_names)]
  } else {
    match.arg(shape)
  }
  display_shape <- unname(shape_names[shape_key])
  acquisition <- match.arg(acquisition)
  noise_type <- match.arg(noise_type)
  n <- as.integer(n)
  k <- as.integer(k)
  samples <- as.integer(samples)
  seed <- as.integer(seed)
  if (is.na(n) || n < 100L || is.na(k) || k < 2L || n < 2L * k ||
      is.na(samples) || samples < 1L || is.na(seed) ||
      !is.finite(noise) || noise < 0 || noise >= 1) {
    stop("Use `n >= 100`, `k >= 2`, `samples >= 1`, and `0 <= noise < 1`.",
         call. = FALSE)
  }

  set.seed(seed)
  candidate_n <- max(8L * n, n + 5000L)
  candidate_xy <- .generate_volume_candidates(candidate_n, acquisition)
  assignment <- .volume_shape_assignment(candidate_xy, display_shape, k)
  truth_id <- assignment$truth
  if (any(tabulate(truth_id, nbins = k) < 2L)) {
    stop("A volumetric class has insufficient candidates in ", display_shape,
         ".", call. = FALSE)
  }

  density <- rep(1, k)
  if (acquisition == "class_imbalanced") {
    density <- sample(exp(seq(log(0.08), log(8), length.out = k)))
  }
  sampling_weight <- density[truth_id]
  minimum_per_class <- min(max(30L, floor(0.003 * n)), floor(n / k))
  mandatory <- unlist(lapply(seq_len(k), function(label) {
    candidates <- which(truth_id == label)
    sample(candidates, min(minimum_per_class, length(candidates)))
  }), use.names = FALSE)
  available <- setdiff(seq_len(candidate_n), mandatory)
  selected <- c(
    mandatory,
    sample(
      available, n - length(mandatory), replace = FALSE,
      prob = sampling_weight[available]
    )
  )
  selected <- sample(selected)

  xy <- candidate_xy[selected, , drop = FALSE]
  truth_id <- truth_id[selected]
  boundary_proximity <- assignment$boundary_proximity[selected]

  sample_id <- integer(n)
  for (label in seq_len(k)) {
    rows <- which(truth_id == label)
    sample_id[rows] <- sample(rep(seq_len(samples), length.out = length(rows)))
  }
  sample_id <- factor(sample_id, levels = seq_len(samples))
  for (sample_index in seq_len(samples)) {
    rows <- which(sample_id == sample_index)
    transform <- .volume_rotation_matrix(
      runif(1L, -0.20, 0.20), runif(1L, -0.12, 0.12)
    )
    xy[rows, ] <- xy[rows, , drop = FALSE] %*% t(transform)
    xy[rows, 1L] <- xy[rows, 1L] + 2.35 * (sample_index - 1L)
    xy[rows, 2L] <- xy[rows, 2L] + runif(1L, -0.05, 0.05)
  }

  labels_id <- truth_id
  corrupted <- logical(n)
  for (sample_index in seq_len(samples)) {
    rows <- which(sample_id == sample_index)
    n_noise <- floor(noise * length(rows))
    if (n_noise < 1L) next
    priority <- if (noise_type == "boundary") {
      order(boundary_proximity[rows])
    } else if (noise_type == "patch") {
      centers <- sample(rows, min(4L, length(rows)))
      nearest <- rep(Inf, length(rows))
      for (center in centers) {
        nearest <- pmin(
          nearest,
          rowSums((sweep(xy[rows, , drop = FALSE], 2L, xy[center, ], "-"))^2)
        )
      }
      order(nearest)
    } else if (noise_type == "region") {
      center <- sample(rows, 1L)
      distance <- rowSums(
        (sweep(xy[rows, , drop = FALSE], 2L, xy[center, ], "-"))^2
      )
      order(distance)
    } else {
      initial <- sample.int(length(rows), n_noise)
      c(initial, setdiff(seq_along(rows), initial))
    }
    corrupt <- rows[.select_class_preserving_corruption(
      priority, n_noise, truth_id[rows], boundary_proximity = boundary_proximity[rows]
    )]

    if (noise_type == "region") {
      center_label <- truth_id[corrupt[1L]]
      target <- .sample_adjacent_volume_label(center_label, assignment$adjacency)
      replacement <- rep.int(target, length(corrupt))
      same <- replacement == truth_id[corrupt]
      if (any(same)) {
        replacement[same] <- vapply(
          truth_id[corrupt[same]], .sample_adjacent_volume_label, integer(1L),
          adjacency = assignment$adjacency
        )
      }
    } else {
      replacement <- vapply(
        truth_id[corrupt], .sample_adjacent_volume_label, integer(1L),
        adjacency = assignment$adjacency
      )
    }
    labels_id[corrupt] <- replacement
    corrupted[corrupt] <- TRUE
  }
  .validate_class_anchors(
    truth_id, labels_id, sample_id, context = "simulate_volumetric_domains()"
  )

  levels <- paste0("cluster", seq_len(k))
  rownames(xy) <- paste0("voxel", seq_len(n))
  colnames(xy) <- c("x", "y", "z")
  truth <- factor(paste0("cluster", truth_id), levels = levels)
  labels <- factor(paste0("cluster", labels_id), levels = levels)
  region <- factor(paste0("region", truth_id), levels = paste0("region", seq_len(k)))
  names(truth) <- names(labels) <- names(sample_id) <- rownames(xy)

  observed_error <- mean(labels != truth)
  expected_error <- sum(vapply(split(seq_len(n), sample_id), function(rows) {
    as.integer(floor(noise * length(rows)))
  }, integer(1L))) / n
  if (any(!is.finite(xy)) || any(table(truth) == 0L) ||
      abs(observed_error - expected_error) > 1 / n) {
    stop("Volumetric simulation validation failed for ", display_shape, ".",
         call. = FALSE)
  }

  region_counts <- table(region)
  boundary <- boundary_proximity <= stats::quantile(
    boundary_proximity, 0.20, names = FALSE
  )
  sparse <- region == names(which.min(region_counts))[1L]
  output <- list(
    xy = xy,
    labels = labels,
    truth = truth,
    region = region,
    samples = sample_id,
    shape = display_shape,
    shape_id = shape_key,
    acquisition = acquisition,
    noise_type = noise_type,
    boundary_proximity = boundary_proximity,
    boundary = boundary,
    sparse = sparse,
    corrupted = corrupted,
    density = stats::setNames(density, levels),
    region_counts = region_counts
  )
  class(output) <- c("spatial_refinement_benchmark", "list")
  output
}

.volume_scale_boundary_distance <- function(distance) {
  reference <- as.numeric(stats::quantile(distance, 0.80, names = FALSE))
  if (!is.finite(reference) || reference <= 0) reference <- max(distance, 1e-8)
  pmin(0.5, 0.5 * distance / max(reference, 1e-8))
}

.volume_ordered_partition <- function(score, k) {
  cuts <- as.numeric(stats::quantile(
    score, seq(0, 1, length.out = k + 1L), names = FALSE, type = 8
  ))
  cuts <- cummax(cuts + seq_along(cuts) * .Machine$double.eps)
  truth <- pmin(k, pmax(1L, findInterval(score, cuts, all.inside = TRUE)))
  internal <- cuts[seq.int(2L, k)]
  distance <- apply(abs(outer(score, internal, "-")), 1L, min)
  list(
    truth = truth,
    boundary_proximity = .volume_scale_boundary_distance(distance)
  )
}

.volume_ordered_adjacency <- function(k) {
  lapply(seq_len(k), function(label) {
    candidates <- c(label - 1L, label + 1L)
    candidates[candidates >= 1L & candidates <= k]
  })
}

.volume_background_adjacency <- function(k) {
  c(list(seq.int(2L, k)), lapply(seq.int(2L, k), function(label) {
    unique(c(
      1L,
      if (label > 2L) label - 1L else integer(),
      if (label < k) label + 1L else integer()
    ))
  }))
}

.volume_shape_assignment <- function(xy, shape, k = 5L) {
  x <- xy[, 1L]
  y <- xy[, 2L]
  z <- xy[, 3L]
  if (shape == "Concentric shells") {
    radius <- sqrt(x^2 + (y / 0.94)^2 + (z / 0.88)^2)
    score <- radius + 0.055 * sin(4 * atan2(y, x)) * (1 - pmin(radius, 1))
    answer <- .volume_ordered_partition(score, k)
    answer$adjacency <- .volume_ordered_adjacency(k)
    return(answer)
  }
  if (shape == "Warped ellipsoids") {
    warped_x <- x - 0.16 * sin(pi * z) + 0.05 * y * z
    warped_y <- y + 0.12 * cos(pi * x) - 0.04 * x * z
    warped_z <- z - 0.10 * sin(pi * y)
    score <- sqrt(
      (warped_x / 0.62)^2 + (warped_y / 0.82)^2 + (warped_z / 0.50)^2
    )
    answer <- .volume_ordered_partition(score, k)
    answer$adjacency <- .volume_ordered_adjacency(k)
    return(answer)
  }
  if (shape == "Toroidal compartments") {
    radial <- sqrt(x^2 + y^2)
    theta <- atan2(y, x)
    score <- sqrt(
      (radial - 0.52 - 0.05 * sin(3 * theta))^2 +
        (z - 0.07 * cos(2 * theta))^2
    )
    answer <- .volume_ordered_partition(score, k)
    answer$adjacency <- .volume_ordered_adjacency(k)
    return(answer)
  }
  if (shape == "Folded layers") {
    score <- z - 0.27 * sin(pi * x) * cos(pi * y) - 0.10 * x * y +
      0.04 * sin(3 * pi * y)
    answer <- .volume_ordered_partition(score, k)
    answer$adjacency <- .volume_ordered_adjacency(k)
    return(answer)
  }
  if (shape == "Thin folded sheets") {
    folded <- z - 0.20 * sin(pi * x) * cos(pi * y) - 0.06 * x * y
    centers <- seq(-0.50, 0.50, length.out = k - 1L)
    distance <- abs(outer(folded, centers, "-"))
    closest <- max.col(-distance)
    closest_distance <- distance[cbind(seq_along(folded), closest)]
    width <- 0.055
    truth <- rep.int(1L, length(folded))
    inside <- closest_distance < width
    truth[inside] <- closest[inside] + 1L
    return(list(
      truth = truth,
      boundary_proximity = .volume_scale_boundary_distance(
        abs(closest_distance - width)
      ),
      adjacency = .volume_background_adjacency(k)
    ))
  }
  if (shape == "Branching tubes") {
    progress <- (z + 1) / 2
    phases <- seq(0, 2 * pi, length.out = k)[seq_len(k - 1L)]
    centers_x <- vapply(phases, function(phase) {
      0.54 * progress * cos(phase) + 0.07 * sin(2 * pi * z + phase)
    }, numeric(length(z)))
    centers_y <- vapply(phases, function(phase) {
      0.54 * progress * sin(phase) + 0.07 * cos(1.5 * pi * z + phase)
    }, numeric(length(z)))
    distance <- (centers_x - x)^2 + (centers_y - y)^2
    closest <- max.col(-distance)
    closest_distance <- sqrt(distance[cbind(seq_along(x), closest)])
    width <- 0.145
    truth <- rep.int(1L, length(x))
    inside <- closest_distance < width
    truth[inside] <- closest[inside] + 1L
    return(list(
      truth = truth,
      boundary_proximity = .volume_scale_boundary_distance(
        abs(closest_distance - width)
      ),
      adjacency = .volume_background_adjacency(k)
    ))
  }
  if (shape == "Disconnected volumes") {
    centers <- rbind(
      c(-0.58, -0.35, -0.35), c(0.50, 0.42, 0.35),
      c(-0.48, 0.42, 0.30), c(0.55, -0.38, -0.25),
      c(0.00, 0.00, 0.58), c(0.08, -0.02, -0.58),
      c(-0.62, 0.02, -0.05), c(0.62, 0.02, 0.08),
      c(-0.08, 0.62, -0.10), c(0.02, -0.62, 0.12)
    )
    scales <- rbind(
      c(0.22, 0.30, 0.24), c(0.28, 0.20, 0.26),
      c(0.24, 0.27, 0.20), c(0.21, 0.29, 0.25),
      c(0.30, 0.24, 0.18), c(0.27, 0.22, 0.20),
      c(0.20, 0.28, 0.27), c(0.22, 0.25, 0.29),
      c(0.28, 0.19, 0.26), c(0.25, 0.23, 0.22)
    )
    distance <- vapply(seq_len(nrow(centers)), function(index) {
      ((x - centers[index, 1L]) / scales[index, 1L])^2 +
        ((y - centers[index, 2L]) / scales[index, 2L])^2 +
        ((z - centers[index, 3L]) / scales[index, 3L])^2
    }, numeric(length(x)))
    closest <- max.col(-distance)
    ordered <- t(apply(distance, 1L, sort, partial = 2L))
    gap <- ordered[, 2L] - ordered[, 1L]
    return(list(
      truth = ((closest - 1L) %% k) + 1L,
      boundary_proximity = .volume_scale_boundary_distance(gap),
      adjacency = lapply(seq_len(k), function(label) setdiff(seq_len(k), label))
    ))
  }
  if (shape == "Helical channels") {
    phases <- seq(0, 2 * pi, length.out = k)[seq_len(k - 1L)]
    centers_x <- vapply(phases, function(phase) {
      0.44 * cos(2.4 * pi * z + phase) * (0.85 - 0.10 * z^2)
    }, numeric(length(z)))
    centers_y <- vapply(phases, function(phase) {
      0.44 * sin(2.4 * pi * z + phase) * (0.85 - 0.10 * z^2)
    }, numeric(length(z)))
    distance <- (centers_x - x)^2 + (centers_y - y)^2
    closest <- max.col(-distance)
    closest_distance <- sqrt(distance[cbind(seq_along(x), closest)])
    width <- 0.125
    truth <- rep.int(1L, length(x))
    inside <- closest_distance < width
    truth[inside] <- closest[inside] + 1L
    return(list(
      truth = truth,
      boundary_proximity = .volume_scale_boundary_distance(
        abs(closest_distance - width)
      ),
      adjacency = .volume_background_adjacency(k)
    ))
  }
  stop("Unknown volumetric shape: ", shape, call. = FALSE)
}

.generate_volume_candidates <- function(n, acquisition) {
  blocks <- list()
  total <- 0L
  block <- 0L
  plane_centers <- c(-0.90, -0.72, -0.49, -0.23, 0.04, 0.31, 0.59, 0.86)
  plane_probability <- c(0.06, 0.13, 0.18, 0.08, 0.21, 0.15, 0.12, 0.07)
  while (total < n) {
    block <- block + 1L
    draw <- max(4000L, as.integer(ceiling(2.8 * (n - total))))
    x <- runif(draw, -1, 1)
    y <- runif(draw, -1, 1)
    z <- if (acquisition == "irregular_z") {
      sample(plane_centers, draw, replace = TRUE, prob = plane_probability) +
        rnorm(draw, sd = 0.012)
    } else {
      runif(draw, -1, 1)
    }
    keep <- x^2 + (y / 0.96)^2 + (z / 0.92)^2 <= 1
    blocks[[block]] <- cbind(x[keep], y[keep], z[keep])
    total <- total + sum(keep)
  }
  xy <- do.call(rbind, blocks)
  xy[seq_len(n), , drop = FALSE]
}

.volume_rotation_matrix <- function(angle_z, angle_y) {
  rz <- matrix(c(
    cos(angle_z), -sin(angle_z), 0,
    sin(angle_z), cos(angle_z), 0,
    0, 0, 1
  ), nrow = 3L, byrow = TRUE)
  ry <- matrix(c(
    cos(angle_y), 0, sin(angle_y),
    0, 1, 0,
    -sin(angle_y), 0, cos(angle_y)
  ), nrow = 3L, byrow = TRUE)
  rz %*% ry
}

.sample_adjacent_volume_label <- function(label, adjacency) {
  candidates <- adjacency[[label]]
  if (!length(candidates)) candidates <- setdiff(seq_along(adjacency), label)
  candidates[sample.int(length(candidates), 1L)]
}
