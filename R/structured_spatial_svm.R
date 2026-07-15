#' Refine spatial cluster assignments with marginSVM
#'
#' `refine_spatial_svm()` repairs noisy assignments using overlapping local
#' Nystrom SVMs, graph evidence, and edge-aware total-variation decoding. Tiles
#' and model settings are selected automatically. When `samples` is supplied,
#' each tissue or section is processed independently.
#'
#' @param xy Numeric matrix with two or three spatial coordinates per row.
#' @param labels Initial cluster assignments, one per row of `xy`.
#' @param samples Optional tissue or section identifier, one per row of `xy`.
#' @param backend Computational backend. `"auto"` uses a registered CUDA or
#'   Metal provider when available and otherwise uses the built-in CPU engine.
#' @param workers Number of CPU tile workers. `NULL` selects up to four physical
#'   cores automatically. Accelerator providers use their own execution policy.
#'
#' @return A factor with the same levels and row names as `labels`. Attributes
#'   `confidence`, `margin`, `local_support`, `tiles`, `backend`, `workers`, and
#'   `abstained_samples` contain diagnostics.
#' @export
#' @examples
#' sim <- simulate_gradient_regions(n = 2000, minority = 0.25, samples = 2)
#' refined <- refine_spatial_svm(sim$xy, sim$labels, sim$samples)
#' mean(refined == sim$truth)
refine_spatial_svm <- function(xy, labels, samples = NULL,
                               backend = c("auto", "cpu", "cuda", "metal"),
                               workers = NULL) {
  .refine_spatial_svm_engine(
    xy = xy, labels = labels, samples = samples, backend = backend,
    workers = workers, control = NULL, seed = 1L, verbose = FALSE
  )
}

.marginsvm_defaults <- function(labels, workers, seed) {
  list(
    neighbors = 12L,
    target_tile_size = 5000L,
    overlap = 0.25,
    gamma = 10,
    landmarks = 48L,
    epochs = 8L,
    learning_rate = 0.01,
    lambda = 1e-4,
    ramp = 2.5,
    retention = if (nlevels(labels) > 10L) 2.00 else 1.00,
    graph_mix = if (nlevels(labels) > 10L) 0.30 else 0.20,
    probability_floor = if (nlevels(labels) > 10L) 0.005 else 1e-8,
    coherence = if (nlevels(labels) > 10L) 6.0 else 0.0,
    preserve_support = if (nlevels(labels) > 10L) 0.60 else 2.0,
    topology_abstention = if (nlevels(labels) > 10L) 1 else 0,
    adaptive_tiles = 1,
    cross_fitting = 1,
    tv_strength = 0.06,
    tv_iterations = 24L,
    workers = workers,
    seed = seed
  )
}

.refine_spatial_svm_engine <- function(xy, labels, samples = NULL,
                                       backend = c("auto", "cpu", "cuda", "metal"),
                                       workers = NULL, control = NULL, seed = 1L,
                                       verbose = FALSE) {
  backend <- match.arg(backend)
  xy <- as.matrix(xy)
  if (!is.numeric(xy) || !length(xy) || ncol(xy) < 2L || ncol(xy) > 3L ||
      any(!is.finite(xy))) {
    stop("`xy` must be a finite numeric matrix with two or three columns.", call. = FALSE)
  }
  storage.mode(xy) <- "double"
  if (length(labels) != nrow(xy) || anyNA(labels)) {
    stop("`labels` must contain one non-missing assignment per row of `xy`.", call. = FALSE)
  }
  labels <- as.factor(labels)
  if (is.null(samples)) samples <- rep.int(1L, nrow(xy))
  if (length(samples) != nrow(xy) || anyNA(samples)) {
    stop("`samples` must contain one non-missing tissue identifier per row.", call. = FALSE)
  }
  samples <- as.integer(as.factor(samples))

  if (is.null(workers)) {
    cores <- parallel::detectCores(logical = FALSE)
    if (is.na(cores)) cores <- 1L
    workers <- max(1L, min(4L, as.integer(cores) - 1L))
  }
  workers <- as.integer(workers)
  seed <- as.integer(seed)
  if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("`workers` must be `NULL` or one positive integer.", call. = FALSE)
  }
  if (length(seed) != 1L || is.na(seed)) {
    stop("Internal `seed` must be one finite integer.", call. = FALSE)
  }

  settings <- .marginsvm_defaults(labels, workers, seed)
  if (is.null(control)) control <- list()
  if (!is.list(control) || (length(control) && is.null(names(control)))) {
    stop("Internal `control` must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(control), names(settings))
  if (length(unknown)) {
    stop("Unknown internal setting: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  settings <- utils::modifyList(settings, control)
  integer_fields <- c("neighbors", "target_tile_size", "landmarks", "epochs",
                      "tv_iterations", "workers", "seed")
  settings[integer_fields] <- lapply(settings[integer_fields], as.integer)
  numeric_fields <- setdiff(names(settings), integer_fields)
  settings[numeric_fields] <- lapply(settings[numeric_fields], as.double)
  valid <- settings$neighbors >= 2L && settings$target_tile_size >= 250L &&
    settings$overlap >= 0 && settings$overlap < 1 && settings$gamma > 0 &&
    settings$landmarks >= 8L && settings$epochs >= 1L &&
    settings$learning_rate > 0 && settings$lambda >= 0 && settings$ramp > 1 &&
    settings$retention >= 0 && settings$graph_mix >= 0 && settings$graph_mix <= 1 &&
    settings$probability_floor > 0 && settings$probability_floor < 1 &&
    settings$coherence >= 0 && settings$preserve_support >= 0 &&
    settings$preserve_support <= 2 && settings$topology_abstention %in% c(0, 1) &&
    settings$adaptive_tiles %in% c(0, 1) && settings$cross_fitting %in% c(0, 1) &&
    settings$tv_strength >= 0 && settings$tv_iterations >= 1L &&
    settings$workers >= 1L &&
    all(vapply(settings, function(value) all(is.finite(value)), logical(1L)))
  if (!valid) stop("Invalid internal marginSVM settings.", call. = FALSE)

  if (nlevels(labels) < 2L) {
    out <- labels
    names(out) <- rownames(xy)
    attr(out, "confidence") <- rep.int(1, nrow(xy))
    attr(out, "margin") <- rep.int(1, nrow(xy))
    attr(out, "local_support") <- rep.int(1, nrow(xy))
    attr(out, "tiles") <- 0L
    attr(out, "abstained_samples") <- integer()
    attr(out, "backend") <- "cpu"
    attr(out, "workers") <- workers
    return(out)
  }

  available <- vapply(c("cuda", "metal"), exists, logical(1L),
                      envir = .spatial_svm_backend_registry, inherits = FALSE)
  if (backend == "auto") {
    preference <- if (Sys.info()[["sysname"]] == "Darwin") {
      c("metal", "cuda", "cpu")
    } else {
      c("cuda", "metal", "cpu")
    }
    backend_used <- preference[preference %in% c("cpu", names(available)[available])][1L]
  } else if (backend == "cpu" || isTRUE(available[[backend]])) {
    backend_used <- backend
  } else {
    warning(backend, " marginSVM backend is unavailable; using CPU.", call. = FALSE)
    backend_used <- "cpu"
  }

  if (backend_used == "cpu") {
    result <- .Call(
      "_marginSVM_structured_spatial_svm_cpp",
      xy, as.integer(labels), samples, settings, isTRUE(verbose),
      PACKAGE = "marginSVM"
    )
  } else {
    result <- .spatial_svm_backend_registry[[backend_used]](
      xy = xy, labels = as.integer(labels), samples = samples, control = settings
    )
    result <- .validate_spatial_svm_backend_result(result, nrow(xy), nlevels(labels))
  }

  refined <- factor(levels(labels)[result$labels], levels = levels(labels))
  names(refined) <- rownames(xy)
  attr(refined, "confidence") <- result$confidence
  attr(refined, "margin") <- result$margin
  attr(refined, "local_support") <- result$local_support
  attr(refined, "tiles") <- result$tiles
  attr(refined, "abstained_samples") <- result$abstained_samples
  attr(refined, "backend") <- backend_used
  attr(refined, "workers") <- settings$workers
  refined
}

.spatial_svm_backend_registry <- new.env(parent = emptyenv())

#' Register a marginSVM accelerator backend
#'
#' This developer interface registers an optional CUDA or Metal provider. A
#' provider receives numeric coordinates, integer labels, integer sample codes,
#' and the frozen settings list, and returns the same diagnostics as the CPU
#' engine.
#'
#' @param backend Either `"cuda"` or `"metal"`.
#' @param provider Backend function, or `NULL` to unregister it.
#' @return Invisibly, the backend name.
#' @export
register_spatial_svm_backend <- function(backend = c("cuda", "metal"), provider) {
  backend <- match.arg(backend)
  if (is.null(provider)) {
    if (exists(backend, envir = .spatial_svm_backend_registry, inherits = FALSE)) {
      rm(list = backend, envir = .spatial_svm_backend_registry)
    }
  } else {
    if (!is.function(provider)) stop("`provider` must be a function or `NULL`.", call. = FALSE)
    assign(backend, provider, envir = .spatial_svm_backend_registry)
  }
  invisible(backend)
}

#' Report marginSVM backend availability
#'
#' @return A data frame reporting availability and implementation source.
#' @export
spatial_svm_backend_capabilities <- function() {
  registered <- unname(vapply(c("cuda", "metal"), exists, logical(1L),
    envir = .spatial_svm_backend_registry, inherits = FALSE))
  data.frame(
    backend = c("cpu", "cuda", "metal"),
    available = c(TRUE, registered),
    source = c("built-in C++", ifelse(registered, "registered provider", "not registered")),
    stringsAsFactors = FALSE
  )
}

.validate_spatial_svm_backend_result <- function(result, n, classes) {
  required <- c("labels", "confidence", "margin", "local_support", "tiles")
  if (!is.list(result) || !all(required %in% names(result))) {
    stop("marginSVM backend returned an incomplete result.", call. = FALSE)
  }
  result$labels <- as.integer(result$labels)
  vectors <- c("labels", "confidence", "margin", "local_support")
  if (any(vapply(result[vectors], length, integer(1L)) != n) ||
      anyNA(result$labels) || any(result$labels < 1L | result$labels > classes)) {
    stop("marginSVM backend returned invalid per-observation values.", call. = FALSE)
  }
  if (is.null(result$abstained_samples)) result$abstained_samples <- integer()
  result
}
