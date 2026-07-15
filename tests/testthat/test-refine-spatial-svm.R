test_that("marginSVM refines a two-dimensional dataset", {
  sim <- simulate_spatial_clusters(n = 700, dimensions = 2, k = 3, noise = 0.15, seed = 11)
  refined <- refine_spatial_svm(sim$xy, sim$labels, workers = 1)

  expect_s3_class(refined, "factor")
  expect_identical(levels(refined), levels(sim$labels))
  expect_length(refined, nrow(sim$xy))
  expect_true(mean(refined == sim$truth) >= mean(sim$labels == sim$truth))
  expect_equal(attr(refined, "backend"), "cpu")
  expect_equal(attr(refined, "workers"), 1L)
})

test_that("three-dimensional tissues are isolated", {
  sim <- simulate_spatial_clusters(
    n = 900, dimensions = 3, k = 4, samples = 3, noise = 0.12, seed = 22
  )
  refined <- refine_spatial_svm(sim$xy, sim$labels, sim$samples, workers = 2)

  expect_length(refined, 900)
  expect_true(attr(refined, "tiles") >= 3L)
  expect_true(all(attr(refined, "confidence") >= 0 &
                  attr(refined, "confidence") <= 1))
  expect_true(all(attr(refined, "margin") >= 0 & attr(refined, "margin") <= 1))
  expect_true(all(attr(refined, "local_support") >= 0 &
                  attr(refined, "local_support") <= 1))
})

test_that("automatic overlapping tiles are used for large inputs", {
  sim <- simulate_gradient_regions(n = 11000, minority = 0.25, samples = 1, seed = 31)
  refined <- refine_spatial_svm(sim$xy, sim$labels, workers = 2)

  expect_gt(attr(refined, "tiles"), 1L)
  expect_gt(mean(refined == sim$truth), mean(sim$labels == sim$truth))
})

test_that("input validation is concise", {
  xy <- matrix(seq_len(20), ncol = 2)
  labels <- rep(c("a", "b"), 5)

  expect_error(refine_spatial_svm(xy[, 1, drop = FALSE], labels), "two or three")
  expect_error(refine_spatial_svm(xy, labels[-1]), "one non-missing")
  expect_error(refine_spatial_svm(xy, labels, samples = labels[-1]), "tissue identifier")
  expect_error(refine_spatial_svm(xy, labels, workers = 0), "positive integer")
})

test_that("one-class input is returned with complete diagnostics", {
  xy <- cbind(x = seq_len(20), y = seq_len(20))
  refined <- refine_spatial_svm(xy, rep("only", 20), workers = 1)

  expect_true(all(refined == "only"))
  expect_equal(attr(refined, "confidence"), rep(1, 20))
  expect_equal(attr(refined, "tiles"), 0L)
  expect_equal(attr(refined, "abstained_samples"), integer())
})

test_that("unavailable accelerators fall back to CPU", {
  sim <- simulate_spatial_clusters(n = 300, seed = 41)
  expect_warning(
    refined <- refine_spatial_svm(sim$xy, sim$labels, backend = "metal", workers = 1),
    "unavailable"
  )
  expect_equal(attr(refined, "backend"), "cpu")
})

test_that("registered accelerator providers follow the production contract", {
  sim <- simulate_spatial_clusters(n = 300, seed = 51)
  provider <- function(xy, labels, samples, control) {
    result <- marginSVM:::.refine_spatial_svm_engine(
      xy, labels, samples, backend = "cpu", workers = 1,
      control = list(topology_abstention = 0), seed = control$seed
    )
    list(
      labels = as.integer(result),
      confidence = attr(result, "confidence"),
      margin = attr(result, "margin"),
      local_support = attr(result, "local_support"),
      tiles = attr(result, "tiles"),
      abstained_samples = attr(result, "abstained_samples")
    )
  }
  register_spatial_svm_backend("metal", provider)
  on.exit(register_spatial_svm_backend("metal", NULL), add = TRUE)

  refined <- refine_spatial_svm(sim$xy, sim$labels, backend = "metal", workers = 1)
  expect_equal(attr(refined, "backend"), "metal")
  expect_true(spatial_svm_backend_capabilities()$available[3L])
})

test_that("the public function exposes only simple execution choices", {
  expect_named(
    formals(refine_spatial_svm),
    c("xy", "labels", "samples", "backend", "workers")
  )
})

test_that("gradient simulator produces the intended A-B-B-C areas", {
  sim <- simulate_gradient_regions(n = 4000, minority = 0.05, samples = 2, seed = 61)
  expect_equal(levels(sim$truth), c("A", "B", "C"))
  expect_equal(length(unique(sim$area)), 4L)
  expect_equal(length(unique(sim$samples)), 2L)
  expect_equal(ncol(sim$xy), 2L)
})
