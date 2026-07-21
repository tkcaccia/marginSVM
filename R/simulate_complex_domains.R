#' Simulate held-out complex spatial domain geometries
#'
#' Generates two- or three-dimensional geometries designed to stress junctions,
#' narrow channels, nonconvex boundaries, and severe regional concentration.
#' These shapes are separate from those in [simulate_spatial_domains()].
#'
#' @param n Number of observations.
#' @param shape Domain geometry.
#' @param density_profile Either `"uniform"` or `"extreme"`.
#' @param noise_type Either `"random"`, `"boundary"`, or `"patch"`.
#' @param noise Fraction of labels to corrupt. Each true class present in a
#'   specimen retains one correct exemplar; a rate that makes this impossible
#'   is rejected.
#' @param samples Number of independent tissues.
#' @param k Number of domains.
#' @param seed Random seed.
#'
#' @return A `spatial_refinement_benchmark`.
#' @export
#' @examples
#' sim <- simulate_complex_spatial_domains(
#'   n = 1000, shape = "tubular_network", density_profile = "extreme"
#' )
#' table(sim$truth)
simulate_complex_spatial_domains <- function(
    n = 50000L,
    shape = c(
      "voronoi_mosaic", "radial_sectors", "checkerboard_junctions",
      "tubular_network", "braided_channels", "shells_3d"
    ),
    density_profile = c("uniform", "extreme"),
    noise_type = c("random", "boundary", "patch"),
    noise = 0.25,
    samples = 1L,
    k = 5L,
    seed = 1L) {
  shape_names <- c(
    voronoi_mosaic = "Voronoi mosaic",
    radial_sectors = "Radial sectors",
    checkerboard_junctions = "Checkerboard junctions",
    tubular_network = "Tubular network",
    braided_channels = "Braided channels",
    shells_3d = "3D shells"
  )
  requested_shape <- shape[1L]
  shape_key <- if (requested_shape %in% unname(shape_names)) {
    names(shape_names)[match(requested_shape, shape_names)]
  } else {
    match.arg(shape)
  }
  display_shape <- unname(shape_names[shape_key])
  density_profile <- match.arg(density_profile)
  noise_type <- match.arg(noise_type)
  n <- as.integer(n)
  samples <- as.integer(samples)
  k <- as.integer(k)
  seed <- as.integer(seed)
  if (is.na(n) || n < 100L || is.na(samples) || samples < 1L ||
      is.na(k) || k < 2L || n < 2L * k || is.na(seed) ||
      !is.finite(noise) || noise < 0 || noise >= 1) {
    stop("Use `n >= 100`, `k >= 2`, `samples >= 1`, and `0 <= noise < 1`.",
         call. = FALSE)
  }

  set.seed(seed)
  candidates <- max(8L * n, n + 2000L)
  if (display_shape == "3D shells") {
    xy <- matrix(rnorm(candidates * 3L), ncol = 3L)
    xy <- xy / sqrt(rowSums(xy^2)) * runif(candidates, 0.12, 1)
  } else {
    xy <- cbind(runif(candidates, -1, 1), runif(candidates, -1, 1))
    inside <- rowSums((sweep(xy, 2L, c(1, 0.88), "/"))^2) <= 1
    xy <- xy[inside, , drop = FALSE]
  }
  truth_id <- .complex_shape_labels(xy, display_shape, k)
  if (any(tabulate(truth_id, nbins = k) == 0L)) {
    stop("A complex-domain class has no candidate observations.", call. = FALSE)
  }
  density <- if (density_profile == "uniform") {
    rep(1, k)
  } else {
    sample(exp(seq(log(0.08), log(8), length.out = k)))
  }
  minimum <- min(max(10L, floor(0.003 * n)), floor(n / k))
  mandatory <- unlist(lapply(seq_len(k), function(label) {
    rows <- which(truth_id == label)
    sample(rows, min(minimum, length(rows)))
  }), use.names = FALSE)
  available <- setdiff(seq_len(nrow(xy)), mandatory)
  selected <- c(
    mandatory,
    sample(
      available, n - length(mandatory),
      prob = density[truth_id[available]]
    )
  )
  selected <- sample(selected)
  xy <- xy[selected, , drop = FALSE]
  truth_id <- truth_id[selected]
  sample_id <- factor(rep(seq_len(samples), length.out = n))
  for (sample in levels(sample_id)) {
    rows <- which(sample_id == sample)
    xy[rows, ] <- sweep(
      xy[rows, , drop = FALSE], 2L, runif(ncol(xy), -0.04, 0.04), "+"
    )
  }

  class_levels <- paste0("cluster", seq_len(k))
  truth <- factor(paste0("cluster", truth_id), levels = class_levels)
  graph_probe <- .potts_icm_labels(
    xy, truth, sample_id,
    control = list(
      iterations = 1L, preserve = 2, consensus = 0.999,
      current_support = 1, margin = 0.999
    )
  )
  boundary_proximity <- attr(graph_probe, "support")
  if (is.null(boundary_proximity)) boundary_proximity <- rep(0.5, n)

  labels_id <- truth_id
  n_noise <- floor(noise * n)
  corrupt <- integer()
  if (n_noise > 0L) {
    priority <- if (noise_type == "boundary") {
      order(boundary_proximity)
    } else if (noise_type == "patch") {
      centers <- sample.int(n, max(2L, ceiling(n_noise / sqrt(n))))
      nearest_patch <- rep(Inf, n)
      for (center in centers) {
        nearest_patch <- pmin(
          nearest_patch,
          rowSums((sweep(xy, 2L, xy[center, ], "-"))^2)
        )
      }
      order(nearest_patch)
    } else {
      initial <- sample.int(n, n_noise)
      c(initial, setdiff(seq_len(n), initial))
    }
    corrupt <- .select_class_preserving_corruption(
      priority, n_noise, truth_id, sample_id, boundary_proximity
    )
    replacement <- sample.int(k, length(corrupt), replace = TRUE)
    same <- replacement == labels_id[corrupt]
    replacement[same] <- replacement[same] %% k + 1L
    labels_id[corrupt] <- replacement
  }
  .validate_class_anchors(
    truth_id, labels_id, sample_id, context = "simulate_complex_spatial_domains()"
  )

  labels <- factor(paste0("cluster", labels_id), levels = class_levels)
  region <- factor(paste0("region", truth_id), levels = paste0("region", seq_len(k)))
  rownames(xy) <- paste0("spot", seq_len(n))
  colnames(xy) <- c("x", "y", "z")[seq_len(ncol(xy))]
  names(labels) <- names(truth) <- names(sample_id) <- rownames(xy)
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
    density_profile = density_profile,
    noise_type = noise_type,
    boundary_proximity = boundary_proximity,
    boundary = boundary,
    sparse = sparse,
    corrupted = labels_id != truth_id,
    density = stats::setNames(density, class_levels),
    region_counts = region_counts,
    density_ratio = max(region_counts) / min(region_counts)
  )
  class(output) <- c("spatial_refinement_benchmark", "list")
  output
}

.complex_shape_labels <- function(xy, shape, k = 5L) {
  x <- xy[, 1L]
  y <- xy[, 2L]
  if (shape == "Voronoi mosaic") {
    centers <- rbind(
      c(-0.55, -0.35), c(0.48, -0.46), c(0.08, 0.02),
      c(-0.42, 0.51), c(0.54, 0.43)
    )
    if (k != 5L) {
      angles <- seq(0, 2 * pi, length.out = k + 1L)[seq_len(k)]
      centers <- cbind(0.58 * cos(angles), 0.50 * sin(angles))
    }
    distance <- vapply(seq_len(k), function(index) {
      (x - centers[index, 1L])^2 + (y - centers[index, 2L])^2
    }, numeric(length(x)))
    return(max.col(-distance))
  }
  if (shape == "Radial sectors") {
    angle <- atan2(y, x) + 0.75 * sqrt(x^2 + y^2) +
      0.18 * sin(4 * atan2(y, x))
    return((floor(((angle + pi) %% (2 * pi)) / (2 * pi) * k) %% k) + 1L)
  }
  if (shape == "Checkerboard junctions") {
    gx <- pmin(5L, floor((x + 1) * 3))
    gy <- pmin(5L, floor((y + 1) * 3))
    return(((gx + 2L * gy) %% k) + 1L)
  }
  if (shape == "Tubular network") {
    if (k != 5L) {
      stop("`tubular_network` currently requires `k = 5`.", call. = FALSE)
    }
    distance <- cbind(
      abs(y - 0.38 * sin(2.5 * pi * x)),
      abs(y + 0.42 * sin(2.1 * pi * x + 0.7)),
      abs(x - 0.34 * sin(2.8 * pi * y + 1.1)),
      abs(x + 0.48 * sin(1.8 * pi * y - 0.4))
    )
    closest <- max.col(-distance)
    answer <- rep(1L, length(x))
    inside <- distance[cbind(seq_along(x), closest)] < 0.075
    answer[inside] <- closest[inside] + 1L
    return(answer)
  }
  if (shape == "Braided channels") {
    score <- y + 0.28 * sin(2.2 * pi * x)
    strand <- floor((score + 1.35) / 2.7 * k)
    crossing <- abs(x) < 0.20
    strand[crossing] <- floor(
      (y[crossing] - 0.25 * sin(8 * pi * x[crossing]) + 1.2) / 2.4 * k
    )
    return(pmax(1L, pmin(k, strand + 1L)))
  }
  radius <- sqrt(rowSums(xy^2)) + 0.10 * sin(3 * atan2(y, x))
  cuts <- stats::quantile(
    radius, seq(0, 1, length.out = k + 1L), names = FALSE
  )
  as.integer(cut(radius, unique(cuts), include.lowest = TRUE, labels = FALSE))
}
