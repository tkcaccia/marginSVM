#' Refine spatial cluster labels with local SVM-like decision boundaries
#'
#' `refine_spatial_svm()` trains local classifiers on spatial coordinates and
#' cluster labels, optionally within spatial tiles and sample strata. The CPU
#' backend is implemented in C++ and approximates an RBF-kernel SVM with random
#' Fourier features plus multiclass hinge-loss optimization. This keeps the
#' method practical for 50k+ observations while retaining nonlinear boundaries.
#'
#' @param xy Numeric matrix of training spatial coordinates. Rows are spots/cells.
#' @param labels Cluster labels for `xy`.
#' @param samples Optional sample/slide identifier for each row of `xy`.
#' @param newdata Optional numeric matrix of coordinates to classify. When `NULL`,
#'   labels are refined for `xy`.
#' @param newsamples Optional sample/slide identifier for each row of `newdata`.
#' @param tiles `"auto"`, `NULL`, or an integer vector with one tile count per
#'   spatial dimension. Automatic tiling is activated only for tissues larger
#'   than `execution$target_tile_size`.
#' @param backend One of `"auto"`, `"cpu"`, `"cuda"`, or `"metal"`. Accelerator
#'   names currently fall back to CPU unless an accelerated build is supplied.
#' @param gamma RBF kernel scale used by the random Fourier features. The
#'   development-benchmark default is 32 after coordinates are normalized
#'   independently within each fitted tissue or tile.
#' @param n_features Number of random Fourier features.
#' @param epochs Number of hinge-loss training passes.
#' @param lambda L2 regularization strength.
#' @param learning_rate Initial learning rate.
#' @param seed Integer seed for reproducible random Fourier features.
#' @param verbose Print sample/tile progress.
#' @param execution Optional named list controlling overlapping tiles. Entries
#'   are `overlap`, `target_tile_size`, `workers`, `blend`, `adaptive`, and
#'   `min_margin`. Defaults combine regular and adaptive tiles with continuous
#'   SVM scores, tapered overlap weights, and conservative margin abstention.
#'
#' @return A factor with the same levels as `labels`.
#' @noRd
#' @examples
#' set.seed(1)
#' sim <- simulate_spatial_clusters(n = 1000, dimensions = 2)
#' refined <- refine_spatial_svm(sim$xy, sim$labels, tiles = c(4, 4))
#' table(refined)
.refine_spatial_svm_legacy <- function(xy,
                               labels,
                               samples = NULL,
                               newdata = NULL,
                               newsamples = NULL,
                               tiles = "auto",
                               backend = c("auto", "cpu", "cuda", "metal"),
                               gamma = NULL,
                               n_features = 256L,
                               epochs = 16L,
                               lambda = 1e-4,
                               learning_rate = 0.005,
                               seed = 1L,
                               verbose = FALSE,
                               execution = NULL) {
  backend <- match.arg(backend)
  xy <- as.matrix(xy)
  storage.mode(xy) <- "double"
  if (!is.numeric(xy) || nrow(xy) == 0L || ncol(xy) < 2L) {
    stop("`xy` must be a non-empty numeric matrix with at least two columns.", call. = FALSE)
  }

  if (length(labels) != nrow(xy)) {
    stop("`labels` must have one value per row of `xy`.", call. = FALSE)
  }
  labels <- as.factor(labels)

  if (is.null(samples)) {
    samples <- factor(rep("sample1", nrow(xy)))
  } else {
    if (length(samples) != nrow(xy)) {
      stop("`samples` must have one value per row of `xy`.", call. = FALSE)
    }
    samples <- as.factor(samples)
  }

  if (is.null(newdata)) {
    newdata_matrix <- matrix(numeric(), nrow = 0L, ncol = ncol(xy))
    newsamples <- factor(character(), levels = levels(samples))
    output_names <- rownames(xy)
  } else {
    newdata_matrix <- as.matrix(newdata)
    storage.mode(newdata_matrix) <- "double"
    if (ncol(newdata_matrix) != ncol(xy)) {
      stop("`newdata` must have the same number of columns as `xy`.", call. = FALSE)
    }
    if (is.null(newsamples)) {
      if (nlevels(samples) > 1L) {
        stop("`newsamples` is required when `samples` has more than one level.", call. = FALSE)
      }
      newsamples <- factor(rep(levels(samples)[1L], nrow(newdata_matrix)), levels = levels(samples))
    } else {
      if (length(newsamples) != nrow(newdata_matrix)) {
        stop("`newsamples` must have one value per row of `newdata`.", call. = FALSE)
      }
      newsamples <- factor(newsamples, levels = levels(samples))
      if (anyNA(newsamples)) {
        stop("`newsamples` contains values not present in `samples`.", call. = FALSE)
      }
    }
    output_names <- rownames(newdata_matrix)
  }

  if (!is.null(tiles) && !identical(tiles, "auto")) {
    tiles <- as.integer(tiles)
    if (length(tiles) != ncol(xy) || anyNA(tiles) || any(tiles < 1L)) {
      stop("`tiles` must be `\"auto\"`, `NULL`, or one positive integer per coordinate dimension.", call. = FALSE)
    }
  }

  if (is.null(gamma)) gamma <- 32
  if (!is.finite(gamma) || gamma <= 0) {
    stop("`gamma` must be a positive finite number.", call. = FALSE)
  }

  execution_defaults <- list(
    overlap = 0.50, target_tile_size = 10000L, workers = NULL,
    blend = "soft", adaptive = "multiscale", min_margin = 0.05
  )
  if (is.null(execution)) execution <- list()
  if (!is.list(execution) || (length(execution) && is.null(names(execution)))) {
    stop("`execution` must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(execution), names(execution_defaults))
  if (length(unknown)) {
    stop("Unknown `execution` setting: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  execution <- utils::modifyList(execution_defaults, execution)
  execution$blend <- match.arg(execution$blend, c("soft", "hard"))
  if (is.logical(execution$adaptive) && length(execution$adaptive) == 1L) {
    execution$adaptive <- if (isTRUE(execution$adaptive)) "adaptive" else "regular"
  } else {
    execution$adaptive <- match.arg(
      execution$adaptive, c("multiscale", "adaptive", "regular")
    )
  }
  if (!is.finite(execution$min_margin) || execution$min_margin < 0 ||
      execution$min_margin >= 1) {
    stop("`execution$min_margin` must be in [0, 1).", call. = FALSE)
  }
  execution$target_tile_size <- as.integer(execution$target_tile_size)
  if (!is.finite(execution$overlap) || execution$overlap < 0 || execution$overlap >= 1 ||
      is.na(execution$target_tile_size) || execution$target_tile_size < 100L) {
    stop("Use `execution$overlap` in [0, 1) and `target_tile_size >= 100`.", call. = FALSE)
  }
  if (!is.null(execution$workers)) {
    execution$workers <- as.integer(execution$workers)
    if (is.na(execution$workers) || execution$workers < 1L) {
      stop("`execution$workers` must be a positive integer.", call. = FALSE)
    }
  }

  if (backend == "metal") {
    warning("Metal SVM backend is unavailable in this build; using CPU.", call. = FALSE)
    backend <- "cpu"
  }
  if (backend == "cuda") {
    warning("CUDA SVM backend is unavailable in this build; using CPU.", call. = FALSE)
    backend <- "cpu"
  }
  backend_code <- switch(backend, auto = 0L, cpu = 1L)

  native_predict <- function(train_xy, train_labels, prediction_xy, block_seed) {
    local_labels <- droplevels(train_labels)
    local_result <- .Call(
      "_SpatialGraphRefine_refine_spatial_svm_cpp",
      train_xy,
      as.integer(local_labels),
      rep.int(1L, nrow(train_xy)),
      prediction_xy,
      rep.int(1L, nrow(prediction_xy)),
      integer(),
      backend_code,
      as.double(gamma),
      as.integer(n_features),
      as.integer(epochs),
      as.double(lambda),
      as.double(learning_rate),
      as.integer(block_seed),
      isTRUE(verbose),
      PACKAGE = "SpatialGraphRefine"
    )
    local_to_global <- match(levels(local_labels), levels(labels))
    global_scores <- matrix(-Inf, nrow = nrow(prediction_xy), ncol = nlevels(labels))
    global_scores[, local_to_global] <- local_result$scores
    list(
      prediction = local_to_global[local_result$labels],
      scores = global_scores
    )
  }

  prediction_xy <- if (is.null(newdata)) xy else newdata_matrix
  prediction_samples <- if (is.null(newdata)) samples else newsamples
  sample_sizes <- table(samples)
  automatic_needs_tiles <- identical(tiles, "auto") && any(sample_sizes > execution$target_tile_size)
  explicit_needs_tiles <- is.integer(tiles) && any(tiles > 1L)
  if (is.null(tiles) || (!automatic_needs_tiles && !explicit_needs_tiles)) {
    result <- .Call(
      "_SpatialGraphRefine_refine_spatial_svm_cpp",
      xy, as.integer(labels), as.integer(samples),
      newdata_matrix, as.integer(newsamples), integer(), backend_code,
      as.double(gamma), as.integer(n_features), as.integer(epochs),
      as.double(lambda), as.double(learning_rate), as.integer(seed),
      isTRUE(verbose), PACKAGE = "SpatialGraphRefine"
    )
    probabilities <- .svm_row_probabilities(result$scores)
    margin <- .svm_probability_margin(probabilities)
    out <- result$labels
    if (is.null(newdata)) {
      retain <- margin < execution$min_margin
      out[retain] <- as.integer(labels)[retain]
    }
    refined <- factor(levels(labels)[out], levels = levels(labels))
    names(refined) <- output_names
    attr(refined, "confidence") <- apply(probabilities, 1L, max)
    attr(refined, "margin") <- margin
    attr(refined, "backend") <- backend
    attr(refined, "workers") <- 1L
    attr(refined, "min_margin") <- execution$min_margin
    return(refined)
  }

  tasks <- .make_svm_tile_tasks(
    xy, labels, samples, prediction_xy, prediction_samples, tiles,
    execution$overlap, execution$target_tile_size, execution$adaptive
  )
  workers <- execution$workers
  if (is.null(workers)) {
    available <- parallel::detectCores(logical = FALSE)
    if (is.na(available)) available <- 1L
    workers <- max(1L, min(length(tasks), available - 1L))
  }
  workers <- max(1L, min(as.integer(workers), length(tasks)))
  if (.Platform$OS.type == "windows" && workers > 1L) {
    warning("Multicore tile execution is unavailable on Windows; using one worker.", call. = FALSE)
    workers <- 1L
  }

  run_task <- function(task) {
    fit <- native_predict(
      xy[task$halo, , drop = FALSE], labels[task$halo],
      prediction_xy[task$prediction, , drop = FALSE],
      seed
    )
    tile_xy <- prediction_xy[task$prediction, , drop = FALSE]
    taper <- sweep(tile_xy, 2L, task$center, "-")
    taper <- abs(sweep(taper, 2L, task$radius, "/"))
    taper <- 1 - taper
    taper[taper < 1e-6] <- 1e-6
    weight <- apply(taper, 1L, prod)^(1 / ncol(tile_xy))
    list(
      core = task$core, prediction_rows = task$prediction,
      prediction = fit$prediction, scores = fit$scores, weight = weight,
      diagnostics = data.frame(
        sample = task$sample, tile = task$tile,
        core_n = length(task$core), prediction_n = length(task$prediction),
        halo_n = length(task$halo),
        scale = if (is.null(task$scale)) "regular" else task$scale,
        depth = if (is.null(task$depth)) NA_integer_ else task$depth,
        impurity = if (is.null(task$impurity)) NA_real_ else task$impurity,
        stringsAsFactors = FALSE
      )
    )
  }
  results <- if (workers > 1L) {
    parallel::mclapply(tasks, run_task, mc.cores = workers, mc.preschedule = TRUE)
  } else {
    lapply(tasks, run_task)
  }
  support <- matrix(0, nrow = nrow(prediction_xy), ncol = nlevels(labels))
  total_weight <- numeric(nrow(prediction_xy))
  diagnostics <- vector("list", length(results))
  for (i in seq_along(results)) {
    result <- results[[i]]
    if (execution$blend == "soft") {
      contribution <- .svm_row_probabilities(result$scores) * result$weight
      weight <- result$weight
    } else {
      contribution <- matrix(0, nrow = length(result$prediction), ncol = nlevels(labels))
      contribution[cbind(seq_along(result$prediction), result$prediction)] <- 1
      weight <- rep.int(1, length(result$prediction))
    }
    support[result$prediction_rows, ] <-
      support[result$prediction_rows, , drop = FALSE] + contribution
    total_weight[result$prediction_rows] <- total_weight[result$prediction_rows] + weight
    diagnostics[[i]] <- results[[i]]$diagnostics
  }
  if (any(total_weight <= 0)) {
    stop("Overlapping SVM tiles did not assign every prediction row.", call. = FALSE)
  }
  probabilities <- support / total_weight
  margin <- .svm_probability_margin(probabilities)
  out <- max.col(probabilities, ties.method = "first")
  if (is.null(newdata)) {
    retain <- margin < execution$min_margin
    out[retain] <- as.integer(labels)[retain]
  }

  refined <- factor(levels(labels)[out], levels = levels(labels))
  names(refined) <- output_names
  attr(refined, "backend") <- backend
  attr(refined, "workers") <- workers
  attr(refined, "tiles") <- do.call(rbind, diagnostics)
  attr(refined, "confidence") <- apply(probabilities, 1L, max)
  attr(refined, "margin") <- margin
  attr(refined, "blend") <- execution$blend
  attr(refined, "min_margin") <- execution$min_margin
  refined
}

.svm_row_probabilities <- function(scores) {
  maximum <- apply(scores, 1L, max)
  probabilities <- exp(scores - maximum)
  probabilities / rowSums(probabilities)
}

.svm_probability_margin <- function(probabilities) {
  if (ncol(probabilities) == 1L) return(rep.int(1, nrow(probabilities)))
  ordered <- t(apply(probabilities, 1L, sort, decreasing = TRUE))
  ordered[, 1L] - ordered[, 2L]
}

.make_svm_tile_tasks <- function(train_xy, train_labels, train_samples,
                                 prediction_xy, prediction_samples,
                                 tiles, overlap, target_tile_size,
                                 adaptive = "multiscale") {
  tasks <- list()
  task_id <- 0L
  for (sample in levels(train_samples)) {
    train_rows <- which(train_samples == sample)
    prediction_rows <- which(prediction_samples == sample)
    if (!length(train_rows) || !length(prediction_rows)) next
    sample_xy <- train_xy[train_rows, , drop = FALSE]
    if (identical(tiles, "auto") && adaptive %in% c("adaptive", "multiscale")) {
      sample_tasks <- .make_adaptive_svm_sample_tasks(
        train_xy, train_labels, train_rows, prediction_xy, prediction_rows,
        sample, overlap, target_tile_size
      )
      for (task in sample_tasks) {
        task_id <- task_id + 1L
        task$id <- task_id
        task$tile <- task_id
        task$scale <- "adaptive"
        tasks[[task_id]] <- task
      }
      if (adaptive == "adaptive") next
    }
    tile_counts <- if (identical(tiles, "auto")) {
      .balanced_spatial_tiles(sample_xy, ceiling(length(train_rows) / target_tile_size))
    } else {
      tiles
    }
    low <- apply(sample_xy, 2L, min)
    high <- apply(sample_xy, 2L, max)
    width <- (high - low) / tile_counts
    width[width <= 0] <- 1
    grid <- expand.grid(lapply(tile_counts, seq_len), KEEP.OUT.ATTRS = FALSE)
    pred_xy <- prediction_xy[prediction_rows, , drop = FALSE]
    for (tile in seq_len(nrow(grid))) {
      position <- as.integer(grid[tile, ])
      tile_low <- low + (position - 1L) * width
      tile_high <- low + position * width
      core <- rep(TRUE, length(prediction_rows))
      prediction_halo <- rep(TRUE, length(prediction_rows))
      halo <- rep(TRUE, length(train_rows))
      for (dimension in seq_len(ncol(train_xy))) {
        if (position[dimension] == 1L) {
          core <- core & pred_xy[, dimension] >= -Inf
        } else {
          core <- core & pred_xy[, dimension] >= tile_low[dimension]
        }
        if (position[dimension] == tile_counts[dimension]) {
          core <- core & pred_xy[, dimension] <= Inf
        } else {
          core <- core & pred_xy[, dimension] < tile_high[dimension]
        }
        halo <- halo & sample_xy[, dimension] >= tile_low[dimension] - overlap * width[dimension]
        halo <- halo & sample_xy[, dimension] <= tile_high[dimension] + overlap * width[dimension]
        prediction_halo <- prediction_halo &
          pred_xy[, dimension] >= tile_low[dimension] - overlap * width[dimension]
        prediction_halo <- prediction_halo &
          pred_xy[, dimension] <= tile_high[dimension] + overlap * width[dimension]
      }
      core_rows <- prediction_rows[core]
      if (!length(core_rows)) next
      prediction_halo <- prediction_halo | core
      task_id <- task_id + 1L
      tasks[[task_id]] <- list(
        id = task_id, sample = sample, tile = tile,
        core = core_rows, prediction = prediction_rows[prediction_halo],
        halo = train_rows[halo], center = (tile_low + tile_high) / 2,
        radius = (.5 + overlap) * width, scale = "regular"
      )
    }
  }
  if (!length(tasks)) stop("SVM tiling produced no prediction tasks.", call. = FALSE)
  tasks
}

.make_adaptive_svm_sample_tasks <- function(train_xy, train_labels, train_rows,
                                            prediction_xy, prediction_rows, sample,
                                            overlap, target_tile_size) {
  dimensions <- ncol(train_xy)
  root_low <- apply(train_xy[train_rows, , drop = FALSE], 2L, min)
  root_high <- apply(train_xy[train_rows, , drop = FALSE], 2L, max)
  minimum_leaf <- max(250L, as.integer(target_tile_size / 4L))
  mandatory_depth <- ceiling(log2(max(1, length(train_rows) / target_tile_size)))
  maximum_depth <- mandatory_depth + 3L

  gini <- function(values) {
    probability <- prop.table(table(values))
    1 - sum(probability^2)
  }
  split_node <- function(rows, low, high, depth) {
    parent_gini <- gini(train_labels[rows])
    best <- NULL
    for (dimension in seq_len(dimensions)) {
      cut <- stats::median(train_xy[rows, dimension])
      left <- rows[train_xy[rows, dimension] <= cut]
      right <- rows[train_xy[rows, dimension] > cut]
      if (length(left) < minimum_leaf || length(right) < minimum_leaf) next
      child_gini <- (length(left) * gini(train_labels[left]) +
        length(right) * gini(train_labels[right])) / length(rows)
      gain <- parent_gini - child_gini
      if (is.null(best) || gain > best$gain) {
        best <- list(dimension = dimension, cut = cut, left = left, right = right, gain = gain)
      }
    }
    must_split <- length(rows) > target_tile_size
    useful_split <- !is.null(best) && best$gain >= 0.025 &&
      length(rows) >= 2L * minimum_leaf
    if (depth >= maximum_depth || (!must_split && !useful_split)) {
      return(list(list(rows = rows, low = low, high = high, depth = depth,
                            impurity = parent_gini)))
    }
    if (is.null(best)) {
      spread <- (high - low) / pmax(root_high - root_low, .Machine$double.eps)
      dimension <- which.max(spread)
      cut <- stats::median(train_xy[rows, dimension])
      best <- list(
        dimension = dimension, cut = cut,
        left = rows[train_xy[rows, dimension] <= cut],
        right = rows[train_xy[rows, dimension] > cut], gain = 0
      )
    }
    if (!length(best$left) || !length(best$right)) {
      return(list(list(rows = rows, low = low, high = high, depth = depth,
                            impurity = parent_gini)))
    }
    left_high <- high
    left_high[best$dimension] <- best$cut
    right_low <- low
    right_low[best$dimension] <- best$cut
    c(
      split_node(best$left, low, left_high, depth + 1L),
      split_node(best$right, right_low, high, depth + 1L)
    )
  }

  leaves <- split_node(train_rows, root_low, root_high, 0L)
  pred_xy <- prediction_xy[prediction_rows, , drop = FALSE]
  owner <- rep.int(NA_integer_, length(prediction_rows))
  for (leaf in seq_along(leaves)) {
    inside <- rep(TRUE, length(prediction_rows))
    for (dimension in seq_len(dimensions)) {
      inside <- inside & pred_xy[, dimension] >= leaves[[leaf]]$low[dimension] &
        pred_xy[, dimension] <= leaves[[leaf]]$high[dimension]
    }
    owner[inside & is.na(owner)] <- leaf
  }
  if (anyNA(owner)) {
    centers <- do.call(rbind, lapply(leaves, function(leaf) (leaf$low + leaf$high) / 2))
    missing <- which(is.na(owner))
    owner[missing] <- apply(pred_xy[missing, , drop = FALSE], 1L, function(point) {
      which.min(rowSums((centers - point)^2))
    })
  }

  tasks <- list()
  for (leaf in seq_along(leaves)) {
    node <- leaves[[leaf]]
    width <- pmax(node$high - node$low, .Machine$double.eps)
    halo_low <- node$low - overlap * width
    halo_high <- node$high + overlap * width
    train_halo <- rep(TRUE, length(train_rows))
    prediction_halo <- rep(TRUE, length(prediction_rows))
    for (dimension in seq_len(dimensions)) {
      train_halo <- train_halo &
        train_xy[train_rows, dimension] >= halo_low[dimension] &
        train_xy[train_rows, dimension] <= halo_high[dimension]
      prediction_halo <- prediction_halo &
        pred_xy[, dimension] >= halo_low[dimension] &
        pred_xy[, dimension] <= halo_high[dimension]
    }
    core <- prediction_rows[owner == leaf]
    if (!length(core)) next
    tasks[[length(tasks) + 1L]] <- list(
      sample = sample, core = core,
      prediction = prediction_rows[prediction_halo | owner == leaf],
      halo = train_rows[train_halo], center = (node$low + node$high) / 2,
      radius = (.5 + overlap) * width, depth = node$depth,
      impurity = node$impurity
    )
  }
  tasks
}

#' Simulate labelled spatial clusters
#'
#' @param n Number of observations.
#' @param dimensions Spatial dimensions, usually 2 or 3.
#' @param k Number of clusters.
#' @param samples Number of independent samples/slides.
#' @param noise Fraction of labels to randomly corrupt.
#' @param seed Random seed.
#'
#' @return A list with `xy`, `labels`, `truth`, and `samples`.
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
  if (dimensions == 3L) {
    centers[, 3L] <- seq(-3, 3, length.out = k)
  }
  truth_id <- sample.int(k, n, replace = TRUE)
  xy <- centers[truth_id, , drop = FALSE] +
    matrix(rnorm(n * dimensions, sd = 0.45), ncol = dimensions)

  for (s in seq_len(samples)) {
    idx <- which(sample_id == s)
    xy[idx, ] <- sweep(xy[idx, , drop = FALSE], 2L, runif(dimensions, -1, 1), "+")
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
