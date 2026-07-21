#' Simulate labelled spatial clusters
#'
#' @param n Number of observations.
#' @param dimensions Spatial dimensions, either 2 or 3.
#' @param k Number of clusters.
#' @param samples Number of independent tissues or sections.
#' @param noise Fraction of labels to corrupt at random. Each true class present
#'   in a specimen retains one correct exemplar; a rate that makes this
#'   impossible is rejected.
#' @param seed Random seed.
#'
#' @return A `spatial_refinement_benchmark` containing coordinates, noisy
#'   labels, truth, sample identifiers, region, boundary, and sparse-region
#'   indicators.
#' @export
#' @examples
#' sim <- simulate_spatial_clusters(n = 500, dimensions = 3, k = 4)
simulate_spatial_clusters <- function(n = 50000L,
                                      dimensions = 2L,
                                      k = 5L,
                                      samples = 1L,
                                      noise = 0.08,
                                      seed = 1L) {
  set.seed(seed)
  dimensions <- as.integer(dimensions)
  k <- as.integer(k)
  samples <- as.integer(samples)
  if (dimensions < 2L || dimensions > 3L) {
    stop("`dimensions` must be 2 or 3.", call. = FALSE)
  }
  if (k < 2L || samples < 1L || n < k || !is.finite(noise) ||
      noise < 0 || noise >= 1) {
    stop("Use `n >= k`, `k >= 2`, `samples >= 1`, and `0 <= noise < 1`.",
         call. = FALSE)
  }

  sample_id <- factor(rep(seq_len(samples), length.out = n))
  angles <- seq(0, 2 * pi, length.out = k + 1L)[seq_len(k)]
  centers <- matrix(0, nrow = k, ncol = dimensions)
  centers[, 1L] <- 4.5 * cos(angles)
  centers[, 2L] <- 4.5 * sin(angles)
  if (dimensions == 3L) centers[, 3L] <- seq(-3, 3, length.out = k)
  truth_id <- sample.int(k, n, replace = TRUE)
  xy <- centers[truth_id, , drop = FALSE] +
    matrix(stats::rnorm(n * dimensions, sd = 0.45), ncol = dimensions)

  offsets <- matrix(0, nrow = samples, ncol = dimensions)
  for (sample in seq_len(samples)) {
    rows <- which(sample_id == sample)
    offsets[sample, ] <- stats::runif(dimensions, -1, 1)
    xy[rows, ] <- sweep(
      xy[rows, , drop = FALSE], 2L, offsets[sample, ], "+"
    )
  }

  boundary_proximity <- numeric(n)
  for (sample in seq_len(samples)) {
    rows <- which(sample_id == sample)
    local_centers <- sweep(centers, 2L, offsets[sample, ], "+")
    distances <- vapply(seq_len(k), function(cluster) {
      rowSums((sweep(xy[rows, , drop = FALSE], 2L,
                     local_centers[cluster, ], "-"))^2)
    }, numeric(length(rows)))
    closest <- max.col(-distances, ties.method = "first")
    closest_distance <- distances[cbind(seq_along(rows), closest)]
    distances[cbind(seq_along(rows), closest)] <- Inf
    second_distance <- apply(distances, 1L, min)
    boundary_proximity[rows] <- sqrt(second_distance) - sqrt(closest_distance)
  }

  labels_id <- truth_id
  n_noise <- floor(noise * n)
  if (n_noise > 0L) {
    initial <- sample.int(n, n_noise)
    corrupt <- .select_class_preserving_corruption(
      c(initial, setdiff(seq_len(n), initial)), n_noise,
      truth_id, sample_id, boundary_proximity
    )
    replacement <- sample.int(k, n_noise, replace = TRUE)
    same <- replacement == truth_id[corrupt]
    replacement[same] <- replacement[same] %% k + 1L
    labels_id[corrupt] <- replacement
  }
  .validate_class_anchors(
    truth_id, labels_id, sample_id, context = "simulate_spatial_clusters()"
  )

  colnames(xy) <- paste0("dim", seq_len(dimensions))
  rownames(xy) <- paste0("spot", seq_len(n))
  class_levels <- paste0("cluster", seq_len(k))
  truth <- factor(paste0("cluster", truth_id), levels = class_levels)
  labels <- factor(paste0("cluster", labels_id), levels = class_levels)
  region <- factor(paste0("region", truth_id), levels = paste0("region", seq_len(k)))
  boundary <- boundary_proximity <= stats::quantile(
    boundary_proximity, 0.20, names = FALSE
  )
  region_count <- table(region)
  sparse <- region == names(which.min(region_count))[1L]
  names(labels) <- names(truth) <- names(sample_id) <- rownames(xy)
  output <- list(
    xy = xy,
    labels = labels,
    truth = truth,
    samples = sample_id,
    region = region,
    boundary_proximity = boundary_proximity,
    boundary = boundary,
    sparse = sparse,
    corrupted = labels != truth,
    region_counts = region_count,
    pattern = "gaussian_clusters",
    density_profile = "uniform"
  )
  class(output) <- c("spatial_refinement_benchmark", "list")
  output
}
