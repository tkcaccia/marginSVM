test_that("refinement returns factor levels and names for training data", {
  sim <- simulate_spatial_clusters(n = 600, dimensions = 2, k = 3, seed = 11)
  refined <- refine_spatial_svm(
    sim$xy,
    sim$labels,
    control = list(target_tile_size = 250L, overlap = 0.20, workers = 2,
                   landmarks = 24L, epochs = 2L, seed = 11L)
  )

  expect_s3_class(refined, "factor")
  expect_identical(levels(refined), levels(sim$labels))
  expect_identical(names(refined), rownames(sim$xy))
  expect_false(anyNA(refined))
  expect_identical(attr(refined, "workers"), 2L)
  expect_gt(attr(refined, "tiles"), 1L)
  expect_true(all(attr(refined, "confidence") >= 0 & attr(refined, "confidence") <= 1))
  expect_true(all(attr(refined, "margin") >= 0 & attr(refined, "margin") <= 1))
})

test_that("multiple tissues are refined independently", {
  sim <- simulate_spatial_clusters(n = 900, dimensions = 3, k = 4, samples = 3, seed = 22)
  refined <- refine_spatial_svm(
    sim$xy,
    sim$labels,
    samples = sim$samples,
    control = list(target_tile_size = 250L, overlap = 0.20, workers = 2,
                   landmarks = 24L, epochs = 2L, seed = 22L)
  )

  expect_length(refined, nrow(sim$xy))
  expect_identical(names(refined), rownames(sim$xy))
  expect_false(anyNA(refined))
  expect_gt(attr(refined, "tiles"), nlevels(sim$samples))
})

test_that("tiled SVM remaps nonconsecutive local factor levels", {
  set.seed(19)
  xy <- rbind(
    cbind(runif(120, 0, 0.4), runif(120)),
    cbind(runif(120, 0.6, 1), runif(120))
  )
  labels <- factor(
    c(rep(c("A", "C"), each = 60), rep(c("B", "C"), each = 60)),
    levels = c("A", "B", "C")
  )
  refined <- refine_spatial_svm(
    xy, labels,
    control = list(target_tile_size = 250L, landmarks = 32L, epochs = 4L,
                   overlap = 0.05, workers = 1L)
  )
  expect_false(anyNA(refined))
  expect_equal(levels(refined), levels(labels))
  expect_false(any(refined[xy[, 1] < 0.4] == "B"))
  expect_false(any(refined[xy[, 1] > 0.6] == "A"))
})

test_that("automatic multiscale tiles cover every tissue", {
  sim <- simulate_gradient_regions(n = 2400, minority = 0.25, samples = 2, seed = 81)
  refined <- refine_spatial_svm(
    sim$xy, sim$labels, sim$samples,
    control = list(target_tile_size = 400L, workers = 1L,
                   landmarks = 32L, epochs = 2L)
  )
  expect_gt(attr(refined, "tiles"), 4L)
  expect_false(anyNA(refined))
})

test_that("structured SVM controls are validated", {
  sim <- simulate_gradient_regions(n = 1200, minority = 0.25, seed = 82)
  expect_error(
    refine_spatial_svm(sim$xy, sim$labels, control = list(unknown = 1)),
    "Unknown.*control"
  )
  expect_error(
    refine_spatial_svm(sim$xy, sim$labels, control = list(graph_mix = 2)),
    "Invalid structured"
  )
})

test_that("structured SVM accelerator providers use the native result contract", {
  sim <- simulate_spatial_clusters(n = 300, dimensions = 2, k = 3, seed = 93)
  expect_warning(
    refine_spatial_svm(sim$xy, sim$labels, backend = "metal"),
    "unavailable"
  )
  provider <- function(xy, labels, samples, control) {
    n <- nrow(xy)
    list(
      labels = labels, confidence = rep.int(1, n), margin = rep.int(1, n),
      local_support = rep.int(1, n), tiles = 1L,
      abstained_samples = unique(samples)
    )
  }
  register_spatial_svm_backend("metal", provider)
  on.exit(register_spatial_svm_backend("metal", NULL), add = TRUE)
  refined <- refine_spatial_svm(sim$xy, sim$labels, backend = "metal")
  expect_identical(as.character(refined), as.character(sim$labels))
  expect_identical(attr(refined, "backend"), "metal")
  expect_true(spatial_svm_backend_capabilities()$available[3L])
})

test_that("gradient regions reproduce the requested mixtures across tissues", {
  sim <- simulate_gradient_regions(n = 8000, minority = 0.05, samples = 2, seed = 71)
  expect_equal(dim(sim$xy), c(8000L, 2L))
  expect_equal(nlevels(sim$samples), 2L)
  area_truth <- vapply(levels(sim$area), function(region) {
    as.character(unique(sim$truth[sim$area == region]))
  }, character(1L))
  expect_identical(unname(area_truth), c("A", "B", "B", "C"))
  observed <- vapply(levels(sim$area), function(region) {
    rows <- sim$area == region
    mean(sim$labels[rows] != sim$truth[rows])
  }, numeric(1L))
  expect_equal(unname(observed), rep(0.05, 4L), tolerance = 0.002)
})

test_that("graph refinement repairs geometric label noise", {
  sim <- simulate_spatial_domains(
    n = 3000,
    pattern = "jagged_stripes",
    k = 5,
    noise = 0.25,
    seed = 44
  )
  refined <- refine_spatial_clusters(sim$xy, sim$labels)
  expect_gt(mean(refined == sim$truth), mean(sim$labels == sim$truth) + 0.08)
  expect_length(attr(refined, "confidence"), nrow(sim$xy))
  expect_identical(attr(refined, "support"), attr(refined, "confidence"))
  expect_true(all(attr(refined, "support") >= 0 & attr(refined, "support") <= 1))
  expect_true(length(attr(refined, "changes")) >= 1L)
})

test_that("graph refinement separates samples and retains names", {
  sim <- simulate_spatial_domains(1200, "rings", samples = 3, seed = 8)
  refined <- refine_spatial_clusters(sim$xy, sim$labels, sim$samples)
  expect_s3_class(refined, "factor")
  expect_identical(names(refined), rownames(sim$xy))
  expect_identical(levels(refined), levels(sim$labels))
  expect_length(attr(refined, "neighbors"), 3L)
  expect_true(all(attr(refined, "neighbors") == 14L))
})

test_that("all geometric simulators return valid labels", {
  patterns <- c(
    "jagged_stripes", "wavy_layers", "rings", "spiral", "branching",
    "lobes", "islands", "disconnected", "thin_layers", "intermixed"
  )
  for (pattern in patterns) {
    sim <- simulate_spatial_domains(500, pattern, k = 4, noise_type = "region", seed = 3)
    expect_equal(nrow(sim$xy), 500)
    expect_false(anyNA(sim$truth))
    expect_equal(nlevels(sim$truth), 4)
    expect_length(sim$boundary_proximity, 500)
    expect_length(sim$corrupted, 500)
  }
  sim3 <- simulate_spatial_domains(500, "layers3d", k = 4, seed = 3)
  expect_equal(ncol(sim3$xy), 3)
})

test_that("clean labels are protected by the discordance gate", {
  sim <- simulate_spatial_domains(3000, "rings", k = 4, noise = 0, seed = 91)
  refined <- refine_spatial_clusters(sim$xy, sim$labels)
  expect_lt(mean(refined != sim$truth), 0.01)
})

test_that("advanced control is validated", {
  sim <- simulate_spatial_domains(300, "lobes", seed = 2)
  expect_error(
    refine_spatial_clusters(sim$xy, sim$labels, control = list(unknown = 1)),
    "Unknown.*control"
  )
  expect_error(
    refine_spatial_clusters(sim$xy, sim$labels, control = list(consensus = 1.5)),
    "Invalid refinement settings"
  )
})

test_that("overlapping 2D tiles preserve ownership and boundary behavior", {
  sim <- simulate_spatial_domains(6000, "wavy_layers", noise = 0.20, seed = 73)
  full <- refine_spatial_clusters(sim$xy, sim$labels)
  tiled <- refine_spatial_clusters(
    sim$xy, sim$labels,
    execution = list(tiles = c(3, 2), overlap = 0.20, workers = 2)
  )
  expect_false(anyNA(tiled))
  expect_identical(names(tiled), rownames(sim$xy))
  expect_gt(mean(tiled == full), 0.99)
  expect_equal(sum(attr(tiled, "tiles")$core_n), nrow(sim$xy))
  expect_true(all(attr(tiled, "tiles")$halo_n >= attr(tiled, "tiles")$core_n))
  expect_identical(attr(tiled, "backend"), "cpu")
  expect_identical(attr(tiled, "workers"), 2L)
})

test_that("overlapping 3D tiles and automatic tiling cover every observation", {
  sim <- simulate_spatial_domains(5000, "layers3d", dimensions = 3, seed = 19)
  tiled <- refine_spatial_clusters(
    sim$xy, sim$labels,
    execution = list(tiles = "auto", target_tile_size = 1000, workers = 1)
  )
  expect_false(anyNA(tiled))
  expect_gt(nrow(attr(tiled, "tiles")), 1L)
  expect_equal(sum(attr(tiled, "tiles")$core_n), nrow(sim$xy))
  expect_equal(ncol(sim$xy), 3L)
})

test_that("accelerator requests are explicit and providers are validated", {
  capabilities <- spatial_backend_capabilities()
  expect_true(capabilities$available[capabilities$backend == "cpu"])
  expect_warning(
    refine_spatial_clusters(
      matrix(runif(600), ncol = 2), factor(rep(1:3, length.out = 300)),
      execution = list(tiles = c(2, 2), backend = "metal", workers = 1)
    ),
    "not registered"
  )
  expect_error(
    refine_spatial_clusters(
      matrix(runif(600), ncol = 2), factor(rep(1:3, length.out = 300)),
      execution = list(tiles = c(2, 2), backend = "cuda", strict_backend = TRUE)
    ),
    "not registered"
  )

  provider <- function(xy, labels, control) {
    SpatialGraphRefine:::.refine_spatial_graph_native(
      xy, labels, rep.int(1L, nrow(xy)), control
    )
  }
  register_spatial_backend("metal", provider)
  on.exit(register_spatial_backend("metal", NULL), add = TRUE)
  sim <- simulate_spatial_domains(1000, "rings", seed = 5)
  accelerated <- refine_spatial_clusters(
    sim$xy, sim$labels,
    execution = list(backend = "metal", workers = 1)
  )
  expect_identical(attr(accelerated, "backend"), "metal")
  expect_false(anyNA(accelerated))
  expect_true(spatial_backend_capabilities()$available[3L])
})

test_that("published direct-refinement rules match source algorithms", {
  set.seed(118)
  xy <- matrix(runif(240), ncol = 2)
  rownames(xy) <- paste0("p", seq_len(nrow(xy)))
  labels <- factor(sample(letters[1:5], nrow(xy), replace = TRUE))
  samples <- factor(rep(1:2, each = nrow(xy) / 2))

  reference <- function(method, k) {
    output <- labels
    for (sample in levels(samples)) {
      rows <- which(samples == sample)
      distances <- as.matrix(dist(xy[rows, , drop = FALSE], upper = TRUE, diag = TRUE))
      for (local in seq_along(rows)) {
        nearest <- order(distances[local, ], seq_along(rows))[-1L][seq_len(k)]
        neighbor_labels <- as.character(labels[rows[nearest]])
        if (method == "graphst") {
          counts <- table(neighbor_labels)
          first_seen <- unique(neighbor_labels)
          output[rows[local]] <- first_seen[which.max(counts[first_seen])]
        } else {
          votes <- c(as.character(labels[rows[local]]), neighbor_labels)
          counts <- table(votes)
          current <- as.character(labels[rows[local]])
          if (counts[current] < k / 2 && max(counts) > k / 2) {
            output[rows[local]] <- names(which.max(counts))
          }
        }
      }
    }
    output
  }

  graphst <- SpatialGraphRefine:::.refine_published_labels(
    xy, labels, samples, method = "graphst", neighbors = 10L
  )
  spagcn <- SpatialGraphRefine:::.refine_published_labels(
    xy, labels, samples, method = "spagcn", neighbors = 6L
  )
  expect_identical(as.character(graphst), as.character(reference("graphst", 10L)))
  expect_identical(as.character(spagcn), as.character(reference("spagcn", 6L)))
})

test_that("experimental marginSVM v2 returns calibrated pointwise diagnostics", {
  sim <- simulate_spatial_domains(
    1800, "thin_layers", k = 5, noise = 0.25,
    noise_type = "boundary", samples = 2, seed = 1901
  )
  refined <- refine_spatial_svm(
    sim$xy, sim$labels, sim$samples,
    control = list(experimental_v2 = 1, workers = 2, seed = 1902)
  )

  expect_s3_class(refined, "factor")
  expect_length(refined, nrow(sim$xy))
  expect_false(anyNA(refined))
  expect_true(all(attr(refined, "trust") >= 0 & attr(refined, "trust") <= 1))
  expect_true(all(attr(refined, "tile_disagreement") >= 0 &
                  attr(refined, "tile_disagreement") <= 1))
  expect_true(all(attr(refined, "perturbation_stability") >= 0 &
                  attr(refined, "perturbation_stability") <= 1))
  expect_true(all(attr(refined, "selective_risk") >= 0 &
                  attr(refined, "selective_risk") <= 1))
  expect_setequal(levels(attr(refined, "decision")),
                  c("retain", "change", "unresolved"))
  expect_length(attr(refined, "protected_component"), nrow(sim$xy))
})

test_that("marginSVM component switches support reproducible ablation", {
  sim <- simulate_spatial_domains(800, "rings", noise = 0.2, seed = 2001)
  regular_in_sample <- refine_spatial_svm(
    sim$xy, sim$labels,
    control = list(adaptive_tiles = 0, cross_fitting = 0, workers = 1, seed = 2002)
  )
  expect_length(regular_in_sample, 800)
  expect_false(anyNA(regular_in_sample))
  expect_identical(attr(regular_in_sample, "control")$adaptive_tiles, 0)
  expect_identical(attr(regular_in_sample, "control")$cross_fitting, 0)
})

test_that("simulated regions support heterogeneous observation concentrations", {
  uniform <- simulate_spatial_domains(
    5000, "jagged_stripes", k = 5, noise = 0,
    density_profile = "uniform", seed = 2101
  )
  strong <- simulate_spatial_domains(
    5000, "jagged_stripes", k = 5, noise = 0,
    density_profile = "strong", seed = 2101
  )
  hotspot <- simulate_spatial_domains(
    5000, "rings", k = 5, noise = 0,
    density_profile = "hotspot", seed = 2102
  )

  expect_equal(nrow(uniform$xy), 5000)
  expect_equal(nrow(strong$xy), 5000)
  expect_equal(nrow(hotspot$xy), 5000)
  expect_identical(strong$density_profile, "strong")
  expect_true(all(strong$region_counts >= 10))
  expect_gt(max(strong$region_counts) / min(strong$region_counts), 4)
  expect_gt(max(hotspot$region_counts) / min(hotspot$region_counts), 2)
  expect_equal(length(strong$region_density), 5)
  expect_error(
    simulate_spatial_domains(500, "rings", k = 5, density_profile = c(1, 2)),
    "one positive finite weight"
  )

  gradient <- simulate_gradient_regions(
    4000, minority = 0.05, samples = 2,
    density_profile = "strong", seed = 2103
  )
  expect_equal(nrow(gradient$xy), 4000)
  expect_equal(as.integer(table(gradient$samples)), c(2000, 2000))
  expect_gt(max(gradient$area_counts) / min(gradient$area_counts), 4)
})
