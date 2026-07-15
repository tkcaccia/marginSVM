#' Refine spatial cluster assignments
#'
#' Repairs locally inconsistent cluster labels while preserving uncertain tissue
#' boundaries. The default settings are inferred from the number of observations,
#' so most analyses only need spatial coordinates and initial cluster labels.
#'
#' @param xy Numeric matrix of 2D or 3D spatial coordinates.
#' @param labels Initial cluster assignments, one per row of `xy`.
#' @param samples Optional slide or sample identifier. Graphs are built separately
#'   for each sample.
#' @param execution Optional named list controlling execution. Supported entries
#'   are `tiles`, `overlap`, `target_tile_size`, `workers`, `backend`, and
#'   `strict_backend`. Use `tiles = "auto"` for overlapping tile-wise execution.
#' @param control Optional named list for advanced use. Supported entries are
#'   `neighbors`, `iterations`, `consensus`, `preserve`, `margin`, and
#'   `weighted`, and `current_support`. These controls are primarily intended for
#'   sensitivity analysis.
#'
#' @return A factor with refined assignments. Attributes `support`, `changes`,
#'   and `neighbors` contain refinement diagnostics. `support` is a local vote
#'   fraction, not a calibrated probability.
#' @keywords internal
#' @noRd
#' @examples
#' sim <- simulate_spatial_domains(n = 2000, pattern = "jagged_stripes")
#' refined <- refine_spatial_clusters(sim$xy, sim$labels)
#' mean(refined == sim$truth)
refine_spatial_clusters <- function(xy, labels, samples = NULL, control = NULL, execution = NULL) {
  xy <- as.matrix(xy)
  storage.mode(xy) <- "double"
  if (!is.numeric(xy) || nrow(xy) < 3L || ncol(xy) < 2L || ncol(xy) > 3L || any(!is.finite(xy))) {
    stop("`xy` must be a finite numeric matrix with 2 or 3 columns and at least 3 rows.", call. = FALSE)
  }
  if (length(labels) != nrow(xy) || anyNA(labels)) {
    stop("`labels` must contain one non-missing value per row of `xy`.", call. = FALSE)
  }
  labels <- as.factor(labels)
  if (nlevels(labels) < 2L) return(labels)

  if (is.null(samples)) {
    samples <- factor(rep("sample1", nrow(xy)))
  } else {
    if (length(samples) != nrow(xy) || anyNA(samples)) {
      stop("`samples` must contain one non-missing value per row of `xy`.", call. = FALSE)
    }
    samples <- as.factor(samples)
  }

  defaults <- list(
    neighbors = NULL,
    iterations = 3L,
    consensus = 0.56,
    preserve = 0.12,
    margin = 0.10,
    weighted = TRUE,
    current_support = NULL
  )
  if (is.null(control)) control <- list()
  if (!is.list(control) || is.null(names(control)) && length(control)) {
    stop("`control` must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown `control` setting: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  settings <- utils::modifyList(defaults, control)
  if (is.null(settings$current_support)) {
    settings$current_support <- 0.25 * min(1, sqrt(12 / nlevels(labels)))
  }
  settings$neighbors <- if (is.null(settings$neighbors)) 0L else as.integer(settings$neighbors)
  settings$iterations <- as.integer(settings$iterations)
  if (is.na(settings$neighbors) || settings$neighbors == 1L || settings$neighbors < 0L ||
      is.na(settings$iterations) || settings$iterations < 1L ||
      !is.finite(settings$consensus) || settings$consensus < 0.5 || settings$consensus >= 1 ||
      !is.finite(settings$preserve) || settings$preserve < 0 ||
      !is.finite(settings$margin) || settings$margin < 0 || settings$margin >= 1 ||
      !is.finite(settings$current_support) || settings$current_support < 0 || settings$current_support > 1 ||
      length(settings$weighted) != 1L || is.na(settings$weighted)) {
    stop("Invalid refinement settings in `control`.", call. = FALSE)
  }

  execution <- .validate_spatial_execution(execution, ncol(xy))
  backend <- .select_spatial_backend(execution$backend, execution$strict_backend)
  tiled <- !is.null(execution$tiles)
  if (!tiled && backend == "cpu") {
    native <- .refine_spatial_graph_native(
      xy, as.integer(labels), as.integer(samples), settings
    )
    return(.format_spatial_refinement(native, labels, rownames(xy), backend, 1L, NULL))
  }
  if (!tiled) execution$tiles <- rep.int(1L, ncol(xy))

  tasks <- .make_spatial_tile_tasks(
    xy, labels, samples, execution$tiles, execution$overlap,
    execution$target_tile_size, settings$neighbors
  )
  workers <- execution$workers
  if (is.null(workers)) {
    available <- parallel::detectCores(logical = FALSE)
    if (is.na(available)) available <- 1L
    workers <- max(1L, min(length(tasks), available - 1L))
  }
  workers <- max(1L, min(as.integer(workers), length(tasks)))
  if (backend != "cpu" && workers > 1L) {
    warning("Accelerator backends process tiles in one R process; setting `workers = 1`.", call. = FALSE)
    workers <- 1L
  }
  if (.Platform$OS.type == "windows" && workers > 1L) {
    warning("Multi-process tile execution currently uses forked workers and is unavailable on Windows; setting `workers = 1`.", call. = FALSE)
    workers <- 1L
  }

  run_task <- function(task) {
    halo <- task$halo
    local_settings <- settings
    local_settings$neighbors <- task$neighbors
    if (backend == "cpu") {
      result <- .refine_spatial_graph_native(
        xy[halo, , drop = FALSE], as.integer(labels)[halo],
        rep.int(1L, length(halo)), local_settings
      )
    } else {
      provider <- .spatial_backend_registry[[backend]]
      result <- provider(
        xy = xy[halo, , drop = FALSE],
        labels = as.integer(labels)[halo],
        control = local_settings
      )
      result <- .validate_backend_result(result, length(halo), nlevels(labels), task$neighbors)
    }
    core_position <- match(task$core, halo)
    list(
      core = task$core,
      labels = as.integer(result)[core_position],
      support = attr(result, "confidence")[core_position],
      changes = sum(as.integer(result)[core_position] != as.integer(labels)[task$core]),
      diagnostics = data.frame(
        sample = as.character(task$sample), tile = task$tile,
        core_n = length(task$core), halo_n = length(halo),
        neighbors = task$neighbors,
        final_changes = sum(as.integer(result)[core_position] != as.integer(labels)[task$core]),
        stringsAsFactors = FALSE
      )
    )
  }
  results <- if (workers > 1L) {
    parallel::mclapply(tasks, run_task, mc.cores = workers, mc.preschedule = TRUE)
  } else {
    lapply(tasks, run_task)
  }
  failed <- vapply(results, inherits, logical(1L), what = "try-error")
  if (any(failed)) {
    stop("Tile worker failed: ", as.character(results[[which(failed)[1L]]]), call. = FALSE)
  }

  output <- as.integer(labels)
  support <- rep(1, nrow(xy))
  max_iterations <- max(vapply(results, function(x) length(x$changes), integer(1L)))
  changes <- integer(max_iterations)
  diagnostics <- vector("list", length(results))
  for (i in seq_along(results)) {
    result <- results[[i]]
    output[result$core] <- result$labels
    support[result$core] <- result$support
    changes[seq_along(result$changes)] <- changes[seq_along(result$changes)] + result$changes
    diagnostics[[i]] <- result$diagnostics
  }
  attr(output, "confidence") <- support
  attr(output, "changes") <- changes
  attr(output, "neighbors") <- vapply(tasks, `[[`, integer(1L), "neighbors")
  tile_diagnostics <- do.call(rbind, diagnostics)
  rownames(tile_diagnostics) <- NULL
  .format_spatial_refinement(output, labels, rownames(xy), backend, workers, tile_diagnostics)
}

.refine_spatial_graph_native <- function(xy, labels, samples, settings) {
  .Call(
    "_marginSVM_refine_spatial_graph_cpp",
    xy,
    labels,
    samples,
    settings$neighbors,
    settings$iterations,
    as.double(settings$consensus),
    as.double(settings$preserve),
    as.double(settings$margin),
    isTRUE(settings$weighted),
    as.double(settings$current_support),
    PACKAGE = "marginSVM"
  )
}

.refine_published_labels <- function(xy, labels, samples, method = c("graphst", "spagcn"), neighbors = NULL) {
  method <- match.arg(method)
  if (is.null(neighbors)) neighbors <- if (method == "graphst") 50L else 6L
  native <- .Call(
    "_marginSVM_direct_refiner_cpp",
    as.matrix(xy), as.integer(labels), as.integer(samples),
    if (method == "graphst") 0L else 1L, as.integer(neighbors),
    PACKAGE = "marginSVM"
  )
  output <- factor(levels(labels)[as.integer(native)], levels = levels(labels))
  names(output) <- rownames(xy)
  attr(output, "neighbors") <- attr(native, "neighbors")
  output
}

.format_spatial_refinement <- function(native, labels, row_names, backend, workers, tiles) {
  diagnostics <- attributes(native)
  refined <- factor(levels(labels)[as.integer(native)], levels = levels(labels))
  names(refined) <- row_names
  attr(refined, "confidence") <- diagnostics$confidence
  attr(refined, "support") <- diagnostics$confidence
  attr(refined, "changes") <- diagnostics$changes
  attr(refined, "neighbors") <- diagnostics$neighbors
  attr(refined, "backend") <- backend
  attr(refined, "workers") <- workers
  attr(refined, "tiles") <- tiles
  refined
}

.validate_spatial_execution <- function(execution, dimensions) {
  defaults <- list(
    tiles = NULL, overlap = 0.15, target_tile_size = 100000L,
    workers = NULL, backend = "auto", strict_backend = FALSE
  )
  if (is.null(execution)) execution <- list()
  if (!is.list(execution) || (length(execution) && is.null(names(execution)))) {
    stop("`execution` must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(execution), names(defaults))
  if (length(unknown)) {
    stop("Unknown `execution` setting: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  out <- utils::modifyList(defaults, execution)
  if (!is.null(out$tiles)) {
    if (identical(out$tiles, "auto")) {
      out$tiles <- "auto"
    } else {
      out$tiles <- as.integer(out$tiles)
      if (length(out$tiles) != dimensions || anyNA(out$tiles) || any(out$tiles < 1L)) {
        stop("`execution$tiles` must be `\"auto\"` or one positive integer per coordinate dimension.", call. = FALSE)
      }
    }
  }
  out$target_tile_size <- as.integer(out$target_tile_size)
  if (!is.finite(out$overlap) || length(out$overlap) != 1L || out$overlap < 0 || out$overlap >= 1 ||
      is.na(out$target_tile_size) || out$target_tile_size < 100L ||
      (!is.null(out$workers) && (length(out$workers) != 1L || is.na(as.integer(out$workers)) || as.integer(out$workers) < 1L)) ||
      length(out$strict_backend) != 1L || is.na(out$strict_backend)) {
    stop("Invalid settings in `execution`.", call. = FALSE)
  }
  out$backend <- match.arg(tolower(out$backend), c("auto", "cpu", "cuda", "metal"))
  out$strict_backend <- isTRUE(out$strict_backend)
  out
}

.balanced_spatial_tiles <- function(xy, target_tiles) {
  dimensions <- ncol(xy)
  counts <- rep.int(1L, dimensions)
  span <- apply(xy, 2L, function(x) diff(range(x)))
  span[span <= 0] <- 1
  while (prod(counts) < target_tiles) {
    dimension <- which.max(span / counts)
    counts[dimension] <- counts[dimension] + 1L
  }
  counts
}

.make_spatial_tile_tasks <- function(xy, labels, samples, tiles, overlap, target_tile_size, requested_neighbors) {
  tasks <- list()
  sample_levels <- levels(samples)
  for (sample in sample_levels) {
    sample_rows <- which(samples == sample)
    sample_xy <- xy[sample_rows, , drop = FALSE]
    sample_tiles <- if (identical(tiles, "auto")) {
      .balanced_spatial_tiles(sample_xy, ceiling(length(sample_rows) / target_tile_size))
    } else {
      tiles
    }
    grid <- expand.grid(lapply(sample_tiles, seq_len), KEEP.OUT.ATTRS = FALSE)
    mins <- apply(sample_xy, 2L, min)
    maxs <- apply(sample_xy, 2L, max)
    widths <- (maxs - mins) / sample_tiles
    widths[widths == 0] <- 1
    sample_classes <- nlevels(droplevels(labels[sample_rows]))
    class_scale <- min(1, 5 / max(1, sample_classes))
    automatic_neighbors <- min(
      31L,
      max(11L, as.integer(round(1.6 * log2(length(sample_rows)) * class_scale)))
    )
    neighbors <- if (requested_neighbors > 0L) requested_neighbors else automatic_neighbors
    for (tile in seq_len(nrow(grid))) {
      position <- as.integer(grid[tile, ])
      low <- mins + (position - 1L) * widths
      high <- mins + position * widths
      core_mask <- rep(TRUE, length(sample_rows))
      halo_mask <- rep(TRUE, length(sample_rows))
      for (dimension in seq_len(ncol(xy))) {
        core_mask <- core_mask & sample_xy[, dimension] >= low[dimension]
        if (position[dimension] < sample_tiles[dimension]) {
          core_mask <- core_mask & sample_xy[, dimension] < high[dimension]
        }
        if (position[dimension] > 1L) {
          halo_mask <- halo_mask & sample_xy[, dimension] >= low[dimension] - overlap * widths[dimension]
        }
        if (position[dimension] < sample_tiles[dimension]) {
          halo_mask <- halo_mask & sample_xy[, dimension] <= high[dimension] + overlap * widths[dimension]
        }
      }
      core <- sample_rows[core_mask]
      halo <- sample_rows[halo_mask]
      if (!length(core)) next
      if (length(halo) < 3L) halo <- core
      tasks[[length(tasks) + 1L]] <- list(
        sample = sample, tile = tile, core = core, halo = halo,
        neighbors = min(neighbors, length(halo) - 1L)
      )
    }
  }
  if (!length(tasks)) stop("Tiling produced no non-empty tile cores.", call. = FALSE)
  tasks
}

.spatial_backend_registry <- new.env(parent = emptyenv())

#' Register an optional accelerator backend
#'
#' Provider packages can register CUDA or Metal implementations without changing
#' the stable user interface. A provider receives `xy`, integer `labels`, and the
#' validated refinement `control` list and returns integer labels with attributes
#' `confidence`, `changes`, and `neighbors`.
#'
#' @param backend Either `"cuda"` or `"metal"`.
#' @param provider Backend function, or `NULL` to unregister it.
#' @return Invisibly, the registered backend name.
#' @keywords internal
#' @noRd
register_spatial_backend <- function(backend = c("cuda", "metal"), provider) {
  backend <- match.arg(backend)
  if (is.null(provider)) {
    if (exists(backend, envir = .spatial_backend_registry, inherits = FALSE)) {
      rm(list = backend, envir = .spatial_backend_registry)
    }
  } else {
    if (!is.function(provider)) stop("`provider` must be a function or `NULL`.", call. = FALSE)
    assign(backend, provider, envir = .spatial_backend_registry)
  }
  invisible(backend)
}

#' Report available spatial refinement backends
#'
#' @return A data frame reporting backend availability and implementation source.
#' @keywords internal
#' @noRd
spatial_backend_capabilities <- function() {
  registered <- unname(vapply(c("cuda", "metal"), exists, logical(1L),
    envir = .spatial_backend_registry, inherits = FALSE))
  data.frame(
    backend = c("cpu", "cuda", "metal"),
    available = c(TRUE, registered),
    source = c("built-in C++", ifelse(registered, "registered provider", "not registered")),
    stringsAsFactors = FALSE
  )
}

.select_spatial_backend <- function(requested, strict) {
  capabilities <- spatial_backend_capabilities()
  available <- capabilities$backend[capabilities$available]
  if (requested == "auto") {
    preference <- if (Sys.info()[["sysname"]] == "Darwin") c("metal", "cuda", "cpu") else c("cuda", "metal", "cpu")
    return(preference[preference %in% available][1L])
  }
  if (requested %in% available) return(requested)
  message <- sprintf("The `%s` backend is not registered; using the built-in CPU backend.", requested)
  if (strict) stop(message, call. = FALSE)
  warning(message, call. = FALSE)
  "cpu"
}

.validate_backend_result <- function(result, n, classes, neighbors) {
  result_attributes <- attributes(result)
  result <- as.integer(result)
  if (length(result) != n || anyNA(result) || any(result < 1L | result > classes)) {
    stop("Accelerator backend returned invalid labels.", call. = FALSE)
  }
  confidence <- result_attributes$confidence
  if (is.null(confidence)) confidence <- rep(NA_real_, n)
  if (length(confidence) != n) stop("Accelerator backend returned invalid support values.", call. = FALSE)
  attr(result, "confidence") <- confidence
  attr(result, "changes") <- if (is.null(result_attributes$changes)) integer() else result_attributes$changes
  attr(result, "neighbors") <- if (is.null(result_attributes$neighbors)) neighbors else result_attributes$neighbors
  result
}

#' Simulate spatial tissue domains with challenging geometry
#'
#' Generates irregular 2D tissue sections or 3D tissue volumes with known domain
#' labels and corrupted initial cluster assignments.
#'
#' @param n Number of observations.
#' @param pattern Domain geometry. See Details.
#' @param k Number of tissue domains.
#' @param noise Fraction of initial labels to corrupt.
#' @param dimensions Either 2 or 3. Pattern `layers3d` requires 3.
#' @param samples Number of independent slides or samples.
#' @param noise_type One of `"random"`, `"boundary"`, `"patch"`, or `"region"`.
#' @param feature_scale Relative width of islands and thin layers. Values below
#'   one create more difficult sub-neighborhood structures.
#' @param density_profile Relative observation concentration across tissue
#'   regions. Use `"uniform"`, `"moderate"`, `"strong"`, `"extreme"`,
#'   `"hotspot"`, or a positive numeric vector of length `k`. Character
#'   profiles permute region weights reproducibly to avoid tying density to a
#'   particular class identifier.
#' @param seed Random seed.
#'
#' @details Available geometries include jagged and wavy layers, concentric
#' rings, spiral arms, branching sectors, lobes, rare islands, disconnected
#' domains, thin layers, interleaved microdomains, and curved 3D layers.
#'
#' @return A list with `xy`, noisy `labels`, `truth`, `samples`, and `pattern`.
#' @export
simulate_spatial_domains <- function(n = 50000L,
                                     pattern = c(
                                       "jagged_stripes", "wavy_layers", "rings",
                                       "spiral", "branching", "lobes", "islands",
                                       "disconnected", "thin_layers", "intermixed",
                                       "layers3d"
                                     ),
                                     k = 5L,
                                     noise = 0.20,
                                     dimensions = if (identical(pattern[1L], "layers3d")) 3L else 2L,
                                     samples = 1L,
                                     noise_type = c("random", "boundary", "patch", "region"),
                                     feature_scale = 1,
                                     density_profile = c("uniform", "moderate", "strong",
                                                         "extreme", "hotspot"),
                                     seed = 1L) {
  pattern <- match.arg(pattern)
  noise_type <- match.arg(noise_type)
  requested_n <- as.integer(n)
  k <- as.integer(k)
  dimensions <- as.integer(dimensions)
  samples <- as.integer(samples)
  if (requested_n < 100L || k < 2L || samples < 1L || dimensions < 2L || dimensions > 3L ||
      !is.finite(noise) || noise < 0 || noise >= 1 ||
      !is.finite(feature_scale) || feature_scale <= 0) {
    stop("Use `n >= 100`, `k >= 2`, 2 or 3 dimensions, and `0 <= noise < 1`.", call. = FALSE)
  }
  if (pattern == "layers3d" && dimensions != 3L) {
    stop("`layers3d` requires `dimensions = 3`.", call. = FALSE)
  }
  if (dimensions == 3L && pattern != "layers3d") {
    stop("Only `layers3d` currently supports three dimensions.", call. = FALSE)
  }

  density_name <- if (is.character(density_profile)) {
    match.arg(density_profile)
  } else {
    "custom"
  }
  if (density_name == "custom") {
    if (!is.numeric(density_profile) || length(density_profile) != k ||
        any(!is.finite(density_profile)) || any(density_profile <= 0)) {
      stop("Numeric `density_profile` must contain one positive finite weight per region.",
           call. = FALSE)
    }
    region_density <- as.numeric(density_profile)
  } else {
    region_density <- switch(
      density_name,
      uniform = rep(1, k),
      moderate = exp(seq(log(0.55), log(1.8), length.out = k)),
      strong = exp(seq(log(0.25), log(4), length.out = k)),
      extreme = exp(seq(log(0.08), log(8), length.out = k)),
      hotspot = exp(seq(log(0.40), log(2.5), length.out = k))
    )
  }

  set.seed(seed)
  if (density_name != "uniform" && density_name != "custom") {
    region_density <- sample(region_density)
  }
  n <- if (density_name == "uniform" && all(region_density == region_density[1L])) {
    requested_n
  } else {
    max(requested_n + 1000L, 6L * requested_n)
  }
  theta <- runif(n, -pi, pi)
  edge <- 1 + 0.14 * sin(3 * theta) + 0.07 * cos(5 * theta)
  radius <- sqrt(runif(n)) * edge
  x <- radius * cos(theta)
  y <- 0.82 * radius * sin(theta)
  boundary_score <- NULL
  boundary_proximity <- rep(0.5, n)
  score_cuts <- NULL

  make_quantile_labels <- function(score) {
    cuts <- unique(stats::quantile(score, seq(0, 1, length.out = k + 1L), names = FALSE))
    score_cuts <<- cuts
    as.integer(cut(score, cuts, include.lowest = TRUE, labels = FALSE))
  }

  score_boundary_proximity <- function(score) {
    internal <- score_cuts[-c(1L, length(score_cuts))]
    if (!length(internal)) return(rep(0.5, length(score)))
    scale <- max(diff(range(score)) / k, .Machine$double.eps)
    pmin(0.5, apply(abs(outer(score, internal, "-")), 1L, min) / scale)
  }

  if (pattern == "jagged_stripes") {
    boundary_score <- x + 0.23 * sin(6 * pi * y) + 0.10 * sin(13 * pi * y)
    truth_id <- make_quantile_labels(boundary_score)
    boundary_proximity <- score_boundary_proximity(boundary_score)
  } else if (pattern == "wavy_layers") {
    boundary_score <- y + 0.30 * sin(2.4 * pi * x) - 0.10 * cos(7 * pi * x)
    truth_id <- make_quantile_labels(boundary_score)
    boundary_proximity <- score_boundary_proximity(boundary_score)
  } else if (pattern == "rings") {
    boundary_score <- radius + 0.055 * sin(9 * theta)
    truth_id <- make_quantile_labels(boundary_score)
    boundary_proximity <- score_boundary_proximity(boundary_score)
  } else if (pattern == "spiral") {
    boundary_score <- ((theta + 3.8 * radius + pi) %% (2 * pi)) / (2 * pi)
    truth_id <- pmin(k, floor(boundary_score * k) + 1L)
    phase <- (boundary_score * k) %% 1
    boundary_proximity <- pmin(phase, 1 - phase)
  } else if (pattern == "branching") {
    boundary_score <- ((theta + pi + 0.75 * radius * sin(3 * theta)) %% (2 * pi)) / (2 * pi)
    truth_id <- pmin(k, floor(boundary_score * k) + 1L)
    phase <- (boundary_score * k) %% 1
    boundary_proximity <- pmin(phase, 1 - phase)
  } else if (pattern == "lobes") {
    angles <- seq(-pi, pi, length.out = k + 1L)[seq_len(k)]
    centers <- cbind(0.58 * cos(angles), 0.48 * sin(angles))
    warped <- cbind(x + 0.13 * sin(4 * y), y + 0.10 * sin(3 * x))
    distances <- vapply(seq_len(k), function(i) rowSums((sweep(warped, 2L, centers[i, ], "-"))^2), numeric(n))
    truth_id <- max.col(-distances)
    ordered <- t(apply(distances, 1L, sort, partial = 2L))
    margin <- ordered[, 2L] - ordered[, 1L]
    boundary_proximity <- pmin(0.5, margin / max(stats::quantile(margin, 0.8), .Machine$double.eps) / 2)
  } else if (pattern == "islands") {
    truth_id <- rep(1L, n)
    angles <- seq(-pi, pi, length.out = k)[seq_len(k - 1L)] + 0.3
    centers <- cbind(0.62 * cos(angles), 0.48 * sin(angles))
    radii <- feature_scale * (0.13 + 0.015 * seq_len(k - 1L))
    circle_distance <- matrix(NA_real_, nrow = n, ncol = k - 1L)
    for (i in seq_len(k - 1L)) {
      radial_distance <- sqrt((x - centers[i, 1L])^2 + (y - centers[i, 2L])^2)
      circle_distance[, i] <- abs(radial_distance - radii[i]) / max(radii[i], .Machine$double.eps)
      inside <- radial_distance < radii[i]
      truth_id[inside] <- i + 1L
    }
    boundary_proximity <- pmin(0.5, apply(circle_distance, 1L, min))
  } else if (pattern == "disconnected") {
    n_centers <- 2L * k
    angles <- seq(-pi, pi, length.out = n_centers + 1L)[seq_len(n_centers)]
    centers <- cbind(0.62 * cos(angles), 0.50 * sin(angles))
    distances <- vapply(seq_len(n_centers), function(i) (x - centers[i, 1L])^2 + (y - centers[i, 2L])^2, numeric(n))
    truth_id <- ((max.col(-distances) - 1L) %% k) + 1L
    ordered <- t(apply(distances, 1L, sort, partial = 2L))
    margin <- ordered[, 2L] - ordered[, 1L]
    boundary_proximity <- pmin(0.5, margin / max(stats::quantile(margin, 0.8), .Machine$double.eps) / 2)
  } else if (pattern == "thin_layers") {
    truth_id <- rep(1L, n)
    centers <- seq(-0.55, 0.55, length.out = k - 1L)
    distances <- vapply(centers, function(center) abs(y - center - 0.08 * sin(3 * pi * x)), numeric(n))
    width <- feature_scale * 0.035
    closest <- max.col(-distances)
    inside <- distances[cbind(seq_len(n), closest)] < width
    truth_id[inside] <- closest[inside] + 1L
    boundary_proximity <- pmin(0.5, abs(distances[cbind(seq_len(n), closest)] - width) / max(width, .Machine$double.eps))
  } else if (pattern == "intermixed") {
    gx <- floor((x - min(x)) / diff(range(x)) * 18)
    gy <- floor((y - min(y)) / diff(range(y)) * 18)
    truth_id <- ((gx + 2L * gy) %% k) + 1L
    cell_x <- abs(((x - min(x)) / diff(range(x)) * 18) %% 1 - 0.5)
    cell_y <- abs(((y - min(y)) / diff(range(y)) * 18) %% 1 - 0.5)
    boundary_proximity <- pmin(0.5, 0.5 - pmax(cell_x, cell_y))
  } else {
    z_limit <- sqrt(pmax(0.04, 1 - pmin(radius, 0.98)^2))
    z <- runif(n, -z_limit, z_limit)
    boundary_score <- z + 0.22 * sin(2.5 * pi * x) - 0.17 * cos(2 * pi * y)
    truth_id <- make_quantile_labels(boundary_score)
    boundary_proximity <- score_boundary_proximity(boundary_score)
  }

  xy <- if (dimensions == 3L) cbind(x, y, z) else cbind(x, y)
  if (n != requested_n) {
    sampling_weight <- region_density[truth_id]
    if (density_name == "hotspot") {
      for (region in seq_len(k)) {
        rows <- which(truth_id == region)
        if (!length(rows)) next
        center <- rows[sample.int(length(rows), 1L)]
        squared_distance <- rowSums((sweep(xy[rows, , drop = FALSE], 2L, xy[center, ], "-"))^2)
        local_scale <- max(stats::quantile(squared_distance, 0.35), .Machine$double.eps)
        sampling_weight[rows] <- sampling_weight[rows] *
          (0.08 + exp(-squared_distance / (2 * local_scale)))
      }
    }
    minimum_per_region <- min(
      max(4L, floor(0.002 * requested_n)),
      floor(requested_n / k)
    )
    mandatory <- unlist(lapply(seq_len(k), function(region) {
      candidates <- which(truth_id == region)
      if (!length(candidates)) return(integer())
      sample(candidates, min(minimum_per_region, length(candidates)), replace = FALSE)
    }), use.names = FALSE)
    remaining_candidates <- setdiff(seq_len(n), mandatory)
    selected <- c(
      mandatory,
      sample(remaining_candidates, requested_n - length(mandatory), replace = FALSE,
             prob = sampling_weight[remaining_candidates])
    )
    selected <- sample(selected)
    xy <- xy[selected, , drop = FALSE]
    truth_id <- truth_id[selected]
    boundary_proximity <- boundary_proximity[selected]
    n <- requested_n
  }
  sample_id <- factor(rep(seq_len(samples), length.out = n))
  for (s in seq_len(samples)) {
    idx <- which(sample_id == s)
    xy[idx, ] <- sweep(xy[idx, , drop = FALSE], 2L, runif(dimensions, -0.08, 0.08), "+")
  }

  labels_id <- truth_id
  n_noise <- floor(noise * n)
  if (n_noise > 0L) {
    if (noise_type == "boundary") {
      corrupt <- order(boundary_proximity)[seq_len(n_noise)]
    } else if (noise_type %in% c("patch", "region")) {
      center_rows <- sample.int(n, max(2L, ceiling(n_noise / max(25L, sqrt(n)))))
      patch_distance <- vapply(center_rows, function(i) rowSums((sweep(xy, 2L, xy[i, ], "-"))^2), numeric(n))
      corrupt <- order(apply(patch_distance, 1L, min))[seq_len(n_noise)]
    } else {
      corrupt <- sample.int(n, n_noise)
    }
    replacement <- if (noise_type == "region") {
      rep(sample.int(k, 1L), n_noise)
    } else {
      sample.int(k, n_noise, replace = TRUE)
    }
    same <- replacement == labels_id[corrupt]
    replacement[same] <- replacement[same] %% k + 1L
    labels_id[corrupt] <- replacement
  }

  colnames(xy) <- c("x", "y", if (dimensions == 3L) "z")
  rownames(xy) <- paste0("spot", seq_len(n))
  list(
    xy = xy,
    labels = factor(paste0("cluster", labels_id), levels = paste0("cluster", seq_len(k))),
    truth = factor(paste0("cluster", truth_id), levels = paste0("cluster", seq_len(k))),
    region = factor(paste0("region", truth_id), levels = paste0("region", seq_len(k))),
    samples = sample_id,
    pattern = pattern,
    noise_type = noise_type,
    density_profile = density_name,
    region_density = stats::setNames(region_density, paste0("cluster", seq_len(k))),
    region_counts = table(factor(paste0("cluster", truth_id),
                                 levels = paste0("cluster", seq_len(k)))),
    boundary_proximity = boundary_proximity,
    corrupted = seq_len(n) %in% if (n_noise > 0L) corrupt else integer()
  )
}
