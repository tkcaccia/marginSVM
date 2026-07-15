#' Refine spatial cluster assignments with marginSVM
#'
#' Fits robust nonlinear SVM boundaries in overlapping adaptive tiles and
#' reconciles their probabilities with an edge-aware spatial total-variation
#' model. All neighborhoods, tiling, fitting, blending, and decoding execute in
#' C++; this R function only validates inputs and selects automatic defaults.
#'
#' @param xy Numeric matrix with two or three spatial coordinates per row.
#' @param labels Initial cluster assignments.
#' @param samples Optional tissue or section identifier per row.
#' @param control Optional named list of advanced controls. Most analyses should
#'   use the automatic defaults.
#' @param backend Computational backend. CPU is always available; registered
#'   accelerator providers can be selected when installed.
#' @param verbose Print native tile progress.
#'
#' @return A factor with confidence, margin, local-support, tile, backend, and
#'   worker diagnostics stored as attributes. The experimental v2 path also
#'   returns trust, tile disagreement, perturbation stability, selective risk,
#'   pointwise decision, and protected-component attributes.
#' @export
refine_spatial_svm <- function(xy, labels, samples = NULL, control = NULL,
                               backend = c("auto", "cpu", "cuda", "metal"),
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
  if (nlevels(labels) < 2L) {
    out <- labels
    names(out) <- rownames(xy)
    attr(out, "confidence") <- rep.int(1, nrow(xy))
    attr(out, "margin") <- rep.int(1, nrow(xy))
    attr(out, "local_support") <- rep.int(1, nrow(xy))
    attr(out, "trust") <- rep.int(1, nrow(xy))
    attr(out, "decision") <- factor(rep.int("retain", nrow(xy)),
      levels = c("retain", "change", "unresolved"))
    attr(out, "tiles") <- 0L
    attr(out, "backend") <- "cpu"
    attr(out, "workers") <- 1L
    return(out)
  }
  if (is.null(samples)) samples <- rep.int(1L, nrow(xy))
  if (length(samples) != nrow(xy) || anyNA(samples)) {
    stop("`samples` must contain one non-missing tissue identifier per row.", call. = FALSE)
  }
  samples <- as.integer(as.factor(samples))

  cores <- parallel::detectCores(logical = FALSE)
  if (is.na(cores)) cores <- 1L
  defaults <- list(
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
    experimental_v2 = 0,
    trust_neighbors = 48L,
    anisotropy = 0.25,
    pairwise_specialists = 1,
    pairwise_max = 8L,
    change_threshold = 0.005,
    unresolved_threshold = 0.001,
    tv_strength = 0.06,
    tv_iterations = 24L,
    workers = max(1L, min(4L, as.integer(cores) - 1L)),
    seed = 1L
  )
  if (is.null(control)) control <- list()
  if (!is.list(control) || (length(control) && is.null(names(control)))) {
    stop("`control` must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown `control` setting: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  control <- utils::modifyList(defaults, control)
  integer_fields <- c("neighbors", "target_tile_size", "landmarks", "epochs",
                      "trust_neighbors", "pairwise_max", "tv_iterations", "workers", "seed")
  control[integer_fields] <- lapply(control[integer_fields], as.integer)
  numeric_fields <- setdiff(names(defaults), c(integer_fields))
  control[numeric_fields] <- lapply(control[numeric_fields], as.double)
  valid <- control$neighbors >= 2L && control$target_tile_size >= 250L &&
    control$overlap >= 0 && control$overlap < 1 && control$gamma > 0 &&
    control$landmarks >= 8L && control$epochs >= 1L &&
    control$learning_rate > 0 && control$lambda >= 0 && control$ramp > 1 &&
    control$retention >= 0 && control$graph_mix >= 0 && control$graph_mix <= 1 &&
    control$probability_floor > 0 && control$probability_floor < 1 &&
    control$coherence >= 0 &&
    control$preserve_support >= 0 && control$preserve_support <= 2 &&
    control$topology_abstention %in% c(0, 1) &&
    control$adaptive_tiles %in% c(0, 1) && control$cross_fitting %in% c(0, 1) &&
    control$experimental_v2 %in% c(0, 1) && control$trust_neighbors >= control$neighbors &&
    control$anisotropy >= 0 && control$anisotropy <= 1 &&
    control$pairwise_specialists %in% c(0, 1) && control$pairwise_max >= 0L &&
    control$change_threshold >= 0 && control$change_threshold <= 1 &&
    control$unresolved_threshold >= 0 &&
    control$unresolved_threshold <= control$change_threshold &&
    control$tv_strength >= 0 &&
    control$tv_iterations >= 1L && control$workers >= 1L &&
    all(vapply(control, function(value) all(is.finite(value)), logical(1L)))
  if (!valid) stop("Invalid structured SVM control settings.", call. = FALSE)

  available <- vapply(c("cuda", "metal"), exists, logical(1L),
                      envir = .spatial_svm_backend_registry, inherits = FALSE)
  if (backend == "auto") {
    preference <- if (Sys.info()[["sysname"]] == "Darwin") {
      c("metal", "cuda", "cpu")
    } else c("cuda", "metal", "cpu")
    backend_used <- preference[preference %in% c("cpu", names(available)[available])][1L]
  } else if (backend == "cpu" || isTRUE(available[[backend]])) {
    backend_used <- backend
  } else {
    warning(backend, " structured SVM backend is unavailable; using CPU.", call. = FALSE)
    backend_used <- "cpu"
  }
  if (backend_used == "cpu") {
    result <- .Call(
      "_SpatialGraphRefine_structured_spatial_svm_cpp",
      xy, as.integer(labels), samples, control, isTRUE(verbose),
      PACKAGE = "SpatialGraphRefine"
    )
  } else {
    result <- .spatial_svm_backend_registry[[backend_used]](
      xy = xy, labels = as.integer(labels), samples = samples, control = control
    )
    result <- .validate_spatial_svm_backend_result(result, nrow(xy), nlevels(labels))
  }
  refined <- factor(levels(labels)[result$labels], levels = levels(labels))
  names(refined) <- rownames(xy)
  attr(refined, "confidence") <- result$confidence
  attr(refined, "margin") <- result$margin
  attr(refined, "local_support") <- result$local_support
  if (!is.null(result$trust)) attr(refined, "trust") <- result$trust
  if (!is.null(result$tile_disagreement)) {
    attr(refined, "tile_disagreement") <- result$tile_disagreement
    attr(refined, "perturbation_stability") <- result$perturbation_stability
    attr(refined, "selective_risk") <- result$selective_risk
    attr(refined, "decision") <- factor(
      c("retain", "change", "unresolved")[result$decision + 1L],
      levels = c("retain", "change", "unresolved"))
    attr(refined, "protected_component") <- result$protected_component
  }
  attr(refined, "tiles") <- result$tiles
  attr(refined, "abstained_samples") <- result$abstained_samples
  attr(refined, "backend") <- backend_used
  attr(refined, "workers") <- control$workers
  attr(refined, "control") <- control
  refined
}

.spatial_svm_backend_registry <- new.env(parent = emptyenv())

#' Register a structured SVM accelerator backend
#'
#' A provider receives numeric `xy`, integer `labels`, integer `samples`, and the
#' validated `control` list. It returns a list containing `labels`, `confidence`,
#' `margin`, `local_support`, and `tiles`; `abstained_samples` is optional.
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

#' Report structured SVM backend availability
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
    stop("Structured SVM backend returned an incomplete result.", call. = FALSE)
  }
  result$labels <- as.integer(result$labels)
  vectors <- c("labels", "confidence", "margin", "local_support")
  if (any(vapply(result[vectors], length, integer(1L)) != n) ||
      anyNA(result$labels) || any(result$labels < 1L | result$labels > classes)) {
    stop("Structured SVM backend returned invalid per-observation values.", call. = FALSE)
  }
  if (is.null(result$abstained_samples)) result$abstained_samples <- integer()
  result
}
