#' Construct a spatial refinement benchmark
#'
#' Creates the common data structure used by the simulation, evaluation, and
#' benchmarking functions in fibermargin.
#'
#' @param xy Numeric matrix containing two or three spatial coordinates per
#'   observation.
#' @param labels Initial cluster assignments.
#' @param truth Reference assignments used only for evaluation.
#' @param samples Optional tissue or section identifier.
#' @param boundary Optional logical vector marking boundary observations.
#' @param regions Optional region identifier used to identify the sparsest and
#'   densest regions.
#' @param sparse Optional logical vector marking observations in sparse regions.
#' @param name Optional benchmark name.
#' @param metadata Optional named list with provenance or scenario information.
#'
#' @return An object of class `spatial_refinement_benchmark`.
#' @export
#' @examples
#' sim <- simulate_spatial_domains(n = 500, pattern = "jagged_stripes")
#' bench <- spatial_benchmark(
#'   sim$xy, sim$labels, sim$truth, sim$samples,
#'   boundary = sim$boundary, regions = sim$region
#' )
spatial_benchmark <- function(xy, labels, truth, samples = NULL,
                              boundary = NULL, regions = NULL, sparse = NULL,
                              name = NULL, metadata = list()) {
  xy <- as.matrix(xy)
  if (!is.numeric(xy) || nrow(xy) < 2L || !ncol(xy) %in% c(2L, 3L) ||
      any(!is.finite(xy))) {
    stop("`xy` must be a finite numeric matrix with two or three columns.",
         call. = FALSE)
  }
  storage.mode(xy) <- "double"
  n <- nrow(xy)
  .validate_benchmark_labels(labels, n, "labels")
  .validate_benchmark_labels(truth, n, "truth")

  if (is.null(samples)) samples <- rep.int("sample1", n)
  .validate_benchmark_labels(samples, n, "samples")
  samples <- as.factor(samples)

  boundary <- .validate_benchmark_flag(boundary, n, "boundary")
  sparse <- .validate_benchmark_flag(sparse, n, "sparse")
  if (!is.null(regions)) {
    .validate_benchmark_labels(regions, n, "regions")
    regions <- as.factor(regions)
  }
  if (!is.null(name) && (length(name) != 1L || is.na(name))) {
    stop("`name` must be `NULL` or one non-missing value.", call. = FALSE)
  }
  if (!is.list(metadata)) {
    stop("`metadata` must be a list.", call. = FALSE)
  }

  if (is.null(rownames(xy))) rownames(xy) <- paste0("spot", seq_len(n))
  labels <- as.factor(labels)
  truth <- as.factor(truth)
  names(labels) <- names(truth) <- names(samples) <- rownames(xy)

  structure(
    list(
      xy = xy,
      labels = labels,
      truth = truth,
      samples = samples,
      boundary = boundary,
      regions = regions,
      sparse = sparse,
      name = if (is.null(name)) NULL else as.character(name),
      metadata = metadata
    ),
    class = c("spatial_refinement_benchmark", "list")
  )
}

#' Evaluate a spatial label refinement
#'
#' Computes complementary measures of label recovery, class balance, boundary
#' behavior, sparse-region behavior, and unintended damage. All measures compare
#' a refined assignment with a fixed reference and the same initial assignment.
#'
#' @param truth Reference assignments.
#' @param initial Initial noisy assignments.
#' @param refined Assignments returned by a refinement method.
#' @param boundary Optional logical vector marking boundary observations.
#' @param regions Optional region identifier. The least and most frequent
#'   regions define sparse- and dense-region accuracy when `sparse` is omitted.
#' @param sparse Optional logical vector marking sparse-region observations.
#' @param elapsed Optional elapsed runtime in seconds.
#' @param method Optional method name included in the returned row.
#'
#' @details
#' `correction_recall` is the fraction of initially wrong labels repaired.
#' `damage_rate` is the fraction of initially correct labels made wrong.
#' `changed_precision` is the fraction of changed labels that are correct after
#' refinement. Empty strata are reported as `NA` rather than silently scored as
#' zero. Adjusted Rand index is calculated internally and does not require an
#' additional package.
#'
#' @return A one-row data frame containing the evaluation measures.
#' @export
#' @examples
#' truth <- factor(c("A", "A", "B", "B"))
#' initial <- factor(c("A", "B", "B", "A"))
#' refined <- truth
#' evaluate_spatial_refinement(truth, initial, refined)
evaluate_spatial_refinement <- function(truth, initial, refined,
                                        boundary = NULL, regions = NULL,
                                        sparse = NULL, elapsed = NA_real_,
                                        method = NULL) {
  n <- length(truth)
  .validate_benchmark_labels(truth, n, "truth")
  .validate_benchmark_labels(initial, n, "initial")
  .validate_benchmark_labels(refined, n, "refined")
  if (!n) stop("Evaluation vectors must not be empty.", call. = FALSE)

  truth_value <- as.character(truth)
  initial_value <- as.character(initial)
  refined_value <- as.character(refined)
  boundary <- .validate_benchmark_flag(boundary, n, "boundary")
  sparse <- .validate_benchmark_flag(sparse, n, "sparse")
  if (!is.null(regions)) {
    .validate_benchmark_labels(regions, n, "regions")
    regions <- as.character(regions)
  }
  if (length(elapsed) != 1L || (!is.na(elapsed) &&
      (!is.numeric(elapsed) || !is.finite(elapsed) || elapsed < 0))) {
    stop("`elapsed` must be `NA` or one non-negative finite number.", call. = FALSE)
  }
  if (!is.null(method) && (length(method) != 1L || is.na(method))) {
    stop("`method` must be `NULL` or one non-missing value.", call. = FALSE)
  }

  correct <- refined_value == truth_value
  initial_correct <- initial_value == truth_value
  changed <- refined_value != initial_value
  truth_classes <- unique(truth_value)
  class_recall <- vapply(truth_classes, function(label) {
    mean(correct[truth_value == label])
  }, numeric(1L))

  sparse_mask <- sparse
  dense_mask <- NULL
  if (!is.null(regions)) {
    region_count <- table(regions)
    if (is.null(sparse_mask)) sparse_mask <- regions == names(which.min(region_count))[1L]
    dense_mask <- regions == names(which.max(region_count))[1L]
  }

  decision <- attr(refined, "decision", exact = TRUE)
  unresolved_fraction <- if (is.null(decision)) {
    unresolved <- attr(refined, "unresolved", exact = TRUE)
    if (is.null(unresolved)) 0 else mean(as.logical(unresolved))
  } else {
    mean(as.character(decision) == "unresolved")
  }
  accuracy <- mean(correct)
  initial_accuracy <- mean(initial_correct)
  initial_error <- 1 - initial_accuracy

  output <- data.frame(
    n = n,
    classes = length(truth_classes),
    initial_accuracy = initial_accuracy,
    accuracy = accuracy,
    accuracy_gain = accuracy - initial_accuracy,
    error_reduction = if (initial_error > 0) {
      (initial_error - (1 - accuracy)) / initial_error
    } else {
      NA_real_
    },
    ari = .adjusted_rand_index(truth_value, refined_value),
    macro_recall = mean(class_recall),
    worst_recall = min(class_recall),
    sparse_region_accuracy = .safe_benchmark_mean(correct, sparse_mask),
    dense_region_accuracy = .safe_benchmark_mean(correct, dense_mask),
    boundary_accuracy = .safe_benchmark_mean(correct, boundary),
    interior_accuracy = .safe_benchmark_mean(
      correct, if (is.null(boundary)) NULL else !boundary
    ),
    correction_recall = .safe_benchmark_mean(correct, !initial_correct),
    damage_rate = .safe_benchmark_mean(!correct, initial_correct),
    changed_precision = .safe_benchmark_mean(correct, changed),
    changed_fraction = mean(changed),
    unresolved_fraction = unresolved_fraction,
    seconds = as.double(elapsed),
    stringsAsFactors = FALSE
  )
  if (!is.null(method)) output <- cbind(method = as.character(method), output)
  class(output) <- c("spatial_refinement_metrics", class(output))
  output
}

#' Benchmark spatial label refinement methods
#'
#' Runs one or more refinement functions on identical benchmark inputs, records
#' elapsed time, and evaluates every output with
#' [evaluate_spatial_refinement()].
#'
#' @param data A benchmark object returned by a fibermargin simulator or
#'   [spatial_benchmark()], or a named list of such objects.
#' @param methods Named list of refinement functions. A function must accept
#'   `xy` and `labels`; `samples` is supplied when the function declares it or
#'   accepts `...`.
#' @param include_initial Include the unrefined assignment as method `Initial`.
#' @param seed Integer seed reset before each method run.
#' @param on_error Either `"stop"` or `"record"`. The latter returns an `error`
#'   column and `NA` performance measures for a failed method.
#'
#' @return A data frame with one row per dataset and method.
#' @export
#' @examples
#' sim <- simulate_spatial_domains(n = 500, pattern = "jagged_stripes", seed = 4)
#' methods <- list(identity = function(xy, labels) labels)
#' benchmark_spatial_refiners(sim, methods)
benchmark_spatial_refiners <- function(
    data,
    methods = list(FiberMargin = refine_spatial_labels),
    include_initial = TRUE,
    seed = 1L,
    on_error = c("stop", "record")) {
  on_error <- match.arg(on_error)
  if (!is.list(methods) || !length(methods) || is.null(names(methods)) ||
      any(!nzchar(names(methods))) || anyDuplicated(names(methods)) ||
      !all(vapply(methods, is.function, logical(1L)))) {
    stop("`methods` must be a non-empty named list of functions.", call. = FALSE)
  }
  if ("Initial" %in% names(methods)) {
    stop("`Initial` is reserved for the unrefined assignment.", call. = FALSE)
  }
  if (length(include_initial) != 1L || is.na(include_initial)) {
    stop("`include_initial` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  seed <- as.integer(seed)
  if (length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be one finite integer.", call. = FALSE)
  }

  required <- c("xy", "labels", "truth")
  if (all(required %in% names(data))) {
    datasets <- list(data)
    names(datasets) <- if (!is.null(data$name) && nzchar(data$name)) {
      data$name
    } else {
      "dataset1"
    }
  } else {
    datasets <- data
    if (!is.list(datasets) || !length(datasets)) {
      stop("`data` must contain one or more benchmark objects.", call. = FALSE)
    }
    if (is.null(names(datasets))) names(datasets) <- paste0("dataset", seq_along(datasets))
    empty_names <- !nzchar(names(datasets))
    names(datasets)[empty_names] <- paste0("dataset", which(empty_names))
    names(datasets) <- make.unique(names(datasets))
  }

  rows <- list()
  row_index <- 0L
  for (dataset_index in seq_along(datasets)) {
    dataset <- .normalize_spatial_benchmark(datasets[[dataset_index]])
    dataset_name <- names(datasets)[dataset_index]
    evaluate <- function(prediction, method, elapsed, error = NA_character_) {
      metrics <- evaluate_spatial_refinement(
        truth = dataset$truth,
        initial = dataset$labels,
        refined = prediction,
        boundary = dataset$boundary,
        regions = dataset$regions,
        sparse = dataset$sparse,
        elapsed = elapsed,
        method = method
      )
      metrics <- cbind(dataset = dataset_name, metrics, error = error)
      class(metrics) <- "data.frame"
      metrics
    }

    if (isTRUE(include_initial)) {
      row_index <- row_index + 1L
      rows[[row_index]] <- evaluate(dataset$labels, "Initial", 0)
    }

    for (method_index in seq_along(methods)) {
      method_name <- names(methods)[method_index]
      set.seed(seed + 1000L * dataset_index + method_index)
      started <- proc.time()[["elapsed"]]
      result <- tryCatch(
        .call_spatial_refiner(methods[[method_index]], dataset),
        error = identity
      )
      elapsed <- unname(proc.time()[["elapsed"]] - started)
      row_index <- row_index + 1L
      if (inherits(result, "error")) {
        if (on_error == "stop") {
          stop("Method `", method_name, "` failed on `", dataset_name,
               "`: ", conditionMessage(result), call. = FALSE)
        }
        failed <- evaluate(dataset$labels, method_name, elapsed,
                           conditionMessage(result))
        metric_columns <- setdiff(
          names(failed),
          c("dataset", "method", "n", "classes", "initial_accuracy",
            "seconds", "error")
        )
        failed[metric_columns] <- NA_real_
        rows[[row_index]] <- failed
      } else {
        rows[[row_index]] <- evaluate(result, method_name, elapsed)
      }
    }
  }
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

.normalize_spatial_benchmark <- function(data) {
  if (!is.list(data) || !all(c("xy", "labels", "truth") %in% names(data))) {
    stop("Each benchmark must contain `xy`, `labels`, and `truth`.", call. = FALSE)
  }
  boundary <- data$boundary
  if (is.null(boundary) && !is.null(data$boundary_proximity)) {
    proximity <- data$boundary_proximity
    if (length(proximity) != nrow(data$xy) || any(!is.finite(proximity))) {
      stop("`boundary_proximity` must contain one finite value per row.",
           call. = FALSE)
    }
    boundary <- proximity <= stats::quantile(proximity, 0.20, names = FALSE)
  }
  regions <- data$regions
  if (is.null(regions)) regions <- data$region
  if (is.null(regions)) regions <- data$area
  sparse <- data$sparse
  if (is.null(sparse)) sparse <- data$sparse_region
  if (is.null(sparse)) sparse <- data$sparse_layer

  benchmark <- spatial_benchmark(
    xy = data$xy,
    labels = data$labels,
    truth = data$truth,
    samples = data$samples,
    boundary = boundary,
    regions = regions,
    sparse = sparse,
    name = data$name,
    metadata = if (is.null(data$metadata)) list() else data$metadata
  )
  benchmark
}

.call_spatial_refiner <- function(method, data) {
  declared <- names(formals(method))
  arguments <- list(xy = data$xy, labels = data$labels)
  if ("samples" %in% declared || "..." %in% declared) {
    arguments$samples <- data$samples
  }
  prediction <- do.call(method, arguments)
  .validate_benchmark_labels(prediction, nrow(data$xy), "method output")
  prediction
}

.validate_benchmark_labels <- function(x, n, name) {
  if (length(x) != n || anyNA(x)) {
    stop("`", name, "` must contain one non-missing value per observation.",
         call. = FALSE)
  }
  invisible(x)
}

.validate_benchmark_flag <- function(x, n, name) {
  if (is.null(x)) return(NULL)
  if (length(x) != n || anyNA(x) || !is.logical(x)) {
    stop("`", name, "` must be `NULL` or one non-missing logical value per observation.",
         call. = FALSE)
  }
  x
}

.safe_benchmark_mean <- function(value, mask) {
  if (is.null(mask) || !any(mask)) return(NA_real_)
  mean(value[mask])
}

.class_anchor_indices <- function(truth, samples = NULL,
                                  boundary_proximity = NULL) {
  n <- length(truth)
  if (!n || anyNA(truth)) {
    stop("Internal error: truth labels must be non-empty and non-missing.",
         call. = FALSE)
  }
  if (is.null(samples)) samples <- rep.int("sample1", n)
  if (length(samples) != n || anyNA(samples)) {
    stop("Internal error: samples must match non-missing truth labels.",
         call. = FALSE)
  }
  if (!is.null(boundary_proximity) &&
      (length(boundary_proximity) != n || any(!is.finite(boundary_proximity)))) {
    stop("Internal error: boundary proximity must match the truth labels.",
         call. = FALSE)
  }

  group <- interaction(
    as.character(samples), as.character(truth), drop = TRUE, lex.order = TRUE
  )
  rows_by_group <- split(seq_len(n), group, drop = TRUE)
  unname(vapply(rows_by_group, function(rows) {
    if (is.null(boundary_proximity)) return(rows[[1L]])
    rows[[which.max(boundary_proximity[rows])]]
  }, integer(1L)))
}

.select_class_preserving_corruption <- function(
    priority, count, truth, samples = NULL, boundary_proximity = NULL) {
  n <- length(truth)
  priority <- as.integer(priority)
  count <- as.integer(count)
  if (length(count) != 1L || is.na(count) || count < 0L) {
    stop("Internal error: corruption count must be one non-negative integer.",
         call. = FALSE)
  }
  if (anyNA(priority) || any(priority < 1L | priority > n)) {
    stop("Internal error: corruption priority contains invalid indices.",
         call. = FALSE)
  }
  priority <- priority[!duplicated(priority)]
  if (length(priority) < count) {
    stop(
      "The requested corruption cannot retain one uncorrupted exemplar of every ",
      "true class within each specimen. Reduce the corruption rate.",
      call. = FALSE
    )
  }
  if (!count) return(integer())
  selected <- priority[seq_len(count)]
  remaining <- priority[-seq_len(count)]
  if (is.null(samples)) samples <- rep.int("sample1", n)
  group <- interaction(
    as.character(samples), as.character(truth), drop = TRUE, lex.order = TRUE
  )
  group_size <- table(group)
  anchor <- .class_anchor_indices(truth, samples, boundary_proximity)
  names(anchor) <- as.character(group[anchor])

  repeat {
    selected_size <- table(factor(group[selected], levels = levels(group)))
    exhausted <- names(selected_size)[selected_size == group_size]
    if (!length(exhausted)) break

    for (key in exhausted) {
      protected <- anchor[[key]]
      selected <- selected[selected != protected]
      replacement <- NA_integer_
      for (candidate in remaining) {
        candidate_key <- as.character(group[[candidate]])
        current_size <- sum(as.character(group[selected]) == candidate_key)
        if (current_size + 1L < unname(group_size[[candidate_key]])) {
          replacement <- candidate
          break
        }
      }
      if (is.na(replacement)) {
        stop(
          "The requested corruption cannot retain one uncorrupted exemplar of every ",
          "true class within each specimen. Reduce the corruption rate.",
          call. = FALSE
        )
      }
      selected <- c(selected, replacement)
      remaining <- remaining[remaining != replacement]
    }
  }
  selected
}

.validate_class_anchors <- function(truth, labels, samples = NULL,
                                    context = "simulation") {
  n <- length(truth)
  if (length(labels) != n || anyNA(truth) || anyNA(labels)) {
    stop("Internal error: labels must match non-missing truth labels.", call. = FALSE)
  }
  if (is.null(samples)) samples <- rep.int("sample1", n)
  if (length(samples) != n || anyNA(samples)) {
    stop("Internal error: samples must match non-missing truth labels.", call. = FALSE)
  }
  truth_value <- as.character(truth)
  label_value <- as.character(labels)
  group <- interaction(
    as.character(samples), truth_value, drop = TRUE, lex.order = TRUE
  )
  retained <- vapply(split(seq_len(n), group, drop = TRUE), function(rows) {
    any(label_value[rows] == truth_value[rows])
  }, logical(1L))
  if (!all(retained)) {
    stop(context, " failed to retain a correct exemplar for every class and specimen.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.adjusted_rand_index <- function(truth, prediction) {
  contingency <- table(truth, prediction)
  choose_two <- function(x) x * (x - 1) / 2
  observations <- sum(contingency)
  if (observations < 2L) return(1)

  index <- sum(choose_two(contingency))
  row_index <- sum(choose_two(rowSums(contingency)))
  column_index <- sum(choose_two(colSums(contingency)))
  total_pairs <- choose_two(observations)
  expected <- row_index * column_index / total_pairs
  maximum <- 0.5 * (row_index + column_index)
  denominator <- maximum - expected
  if (abs(denominator) <= .Machine$double.eps * max(1, abs(maximum))) {
    equivalent <- all(rowSums(contingency > 0) == 1L) &&
      all(colSums(contingency > 0) == 1L)
    return(if (equivalent) 1 else 0)
  }
  (index - expected) / denominator
}
