#' Simulate layered gradient-mixture tissue regions
#'
#' Generates four ordered tissue areas whose reference classes are `A`, `B`,
#' `B`, and `C`. Within the areas, a configurable minority fraction is labeled
#' `B`, `A`, `C`, and `B`, respectively. Independent tissues are generated with
#' separate curved boundaries and coordinate offsets.
#'
#' @param n Total number of observations across all tissues.
#' @param minority Fraction of minority labels in every area.
#' @param dimensions Either 2 or 3 spatial dimensions.
#' @param samples Number of independent tissues.
#' @param seed Random seed.
#' @param curvature Boundary waviness on the unit coordinate scale.
#' @param density_profile Relative observation concentration across the four
#'   tissue areas. Use `"uniform"`, `"moderate"`, `"strong"`, `"extreme"`,
#'   or a positive numeric vector of length four.
#'
#' @return A list containing `xy`, mixed `labels`, reference `truth`, `samples`,
#'   `area`, and normalized `boundary_proximity`.
#' @export
#' @examples
#' sim <- simulate_gradient_regions(n = 4000, minority = 0.05, samples = 2)
#' table(sim$area, sim$labels)
simulate_gradient_regions <- function(n = 50000L,
                                      minority = 0.05,
                                      dimensions = 2L,
                                      samples = 1L,
                                      seed = 1L,
                                      curvature = 0.04,
                                      density_profile = c("uniform", "moderate",
                                                          "strong", "extreme")) {
  requested_n <- as.integer(n)
  dimensions <- as.integer(dimensions)
  samples <- as.integer(samples)
  if (is.na(requested_n) || requested_n < 8L || is.na(samples) || samples < 1L ||
      requested_n < 4L * samples) {
    stop("Use `n >= 8`, `samples >= 1`, and at least four observations per tissue.", call. = FALSE)
  }
  if (!dimensions %in% c(2L, 3L)) {
    stop("`dimensions` must be 2 or 3.", call. = FALSE)
  }
  if (!is.finite(minority) || minority < 0 || minority >= 0.5) {
    stop("`minority` must be in [0, 0.5).", call. = FALSE)
  }
  if (!is.finite(curvature) || curvature < 0 || curvature >= 0.2) {
    stop("`curvature` must be in [0, 0.2).", call. = FALSE)
  }

  density_name <- if (is.character(density_profile)) match.arg(density_profile) else "custom"
  area_density <- if (density_name == "custom") {
    if (!is.numeric(density_profile) || length(density_profile) != 4L ||
        any(!is.finite(density_profile)) || any(density_profile <= 0)) {
      stop("Numeric `density_profile` must contain four positive finite weights.", call. = FALSE)
    }
    as.numeric(density_profile)
  } else {
    switch(density_name,
      uniform = rep(1, 4L),
      moderate = exp(seq(log(0.55), log(1.8), length.out = 4L)),
      strong = exp(seq(log(0.25), log(4), length.out = 4L)),
      extreme = exp(seq(log(0.08), log(8), length.out = 4L)))
  }

  set.seed(seed)
  if (density_name != "uniform" && density_name != "custom") area_density <- sample(area_density)
  n <- if (density_name == "uniform" && all(area_density == area_density[1L])) {
    requested_n
  } else {
    6L * requested_n
  }
  sample_id <- factor(rep(seq_len(samples), length.out = n))
  xy <- matrix(0, nrow = n, ncol = dimensions)
  area <- integer(n)
  boundary_proximity <- numeric(n)

  for (sample in seq_len(samples)) {
    rows <- which(sample_id == sample)
    count <- length(rows)
    x <- runif(count)
    y <- runif(count)
    phase <- runif(1L, 0, 2 * pi)
    score <- y + curvature * sin(2 * pi * x + phase)
    cuts <- as.numeric(stats::quantile(score, c(0, 0.25, 0.5, 0.75, 1), names = FALSE))
    local_area <- pmin(4L, findInterval(score, cuts, all.inside = TRUE))
    local_area[local_area < 1L] <- 1L
    area[rows] <- local_area
    internal <- cuts[2:4]
    boundary_proximity[rows] <- apply(abs(outer(score, internal, "-")), 1L, min)
    xy[rows, 1L] <- x + 1.25 * (sample - 1L)
    xy[rows, 2L] <- y
    if (dimensions == 3L) {
      xy[rows, 3L] <- runif(count) + 0.08 * sin(2 * pi * x)
    }
  }

  if (n != requested_n) {
    requested_sample <- factor(rep(seq_len(samples), length.out = requested_n))
    selected <- integer()
    for (sample in seq_len(samples)) {
      candidates <- which(sample_id == sample)
      take <- sum(requested_sample == sample)
      selected <- c(selected, sample(candidates, take, replace = FALSE,
                                     prob = area_density[area[candidates]]))
    }
    xy <- xy[selected, , drop = FALSE]
    area <- area[selected]
    boundary_proximity <- boundary_proximity[selected]
    sample_id <- factor(sample_id[selected], levels = seq_len(samples))
    n <- requested_n
  }

  truth_names <- c("A", "B", "B", "C")
  minority_names <- c("B", "A", "C", "B")
  truth <- factor(truth_names[area], levels = c("A", "B", "C"))
  labels <- truth
  for (sample in seq_len(samples)) {
    for (region in seq_len(4L)) {
      rows <- which(sample_id == sample & area == region)
      contaminated <- if (minority == 0) integer() else {
        sample(rows, as.integer(round(minority * length(rows))))
      }
      labels[contaminated] <- minority_names[region]
    }
  }

  rownames(xy) <- paste0("spot", seq_len(n))
  colnames(xy) <- c("x", "y", "z")[seq_len(dimensions)]
  names(labels) <- names(truth) <- names(sample_id) <- rownames(xy)
  list(
    xy = xy,
    labels = labels,
    truth = truth,
    samples = sample_id,
    area = factor(area, levels = 1:4, labels = paste0("area", 1:4)),
    density_profile = density_name,
    area_density = stats::setNames(area_density, paste0("area", 1:4)),
    area_counts = table(factor(area, levels = 1:4, labels = paste0("area", 1:4))),
    boundary_proximity = boundary_proximity
  )
}
