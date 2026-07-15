#' Simulate labelled spatial clusters
#'
#' @param n Number of observations.
#' @param dimensions Spatial dimensions, either 2 or 3.
#' @param k Number of clusters.
#' @param samples Number of independent tissues or sections.
#' @param noise Fraction of labels to corrupt at random.
#' @param seed Random seed.
#'
#' @return A list with `xy`, `labels`, `truth`, and `samples`.
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
  if (k < 2L || samples < 1L || n < k) {
    stop("Use `n >= k`, `k >= 2`, and `samples >= 1`.", call. = FALSE)
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

  for (sample in seq_len(samples)) {
    rows <- which(sample_id == sample)
    xy[rows, ] <- sweep(
      xy[rows, , drop = FALSE], 2L, stats::runif(dimensions, -1, 1), "+"
    )
  }

  labels_id <- truth_id
  n_noise <- floor(noise * n)
  if (n_noise > 0L) {
    corrupt <- sample.int(n, n_noise)
    labels_id[corrupt] <- sample.int(k, n_noise, replace = TRUE)
  }

  colnames(xy) <- paste0("dim", seq_len(dimensions))
  rownames(xy) <- paste0("spot", seq_len(n))
  list(
    xy = xy,
    labels = factor(paste0("cluster", labels_id)),
    truth = factor(paste0("cluster", truth_id)),
    samples = sample_id
  )
}
