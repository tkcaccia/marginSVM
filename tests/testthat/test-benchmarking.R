test_that("evaluation measures correction and damage from the same initial labels", {
  truth <- factor(c("A", "A", "B", "B", "C", "C"))
  initial <- factor(c("A", "B", "B", "B", "C", "A"))
  refined <- factor(c("A", "A", "B", "A", "C", "C"))
  boundary <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
  sparse <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)

  metrics <- evaluate_spatial_refinement(
    truth, initial, refined, boundary = boundary, sparse = sparse,
    method = "test"
  )

  expect_identical(metrics$method, "test")
  expect_equal(metrics$initial_accuracy, 4 / 6)
  expect_equal(metrics$accuracy, 5 / 6)
  expect_equal(metrics$accuracy_gain, 1 / 6)
  expect_equal(metrics$correction_recall, 1)
  expect_equal(metrics$damage_rate, 1 / 4)
  expect_equal(metrics$changed_precision, 2 / 3)
  expect_equal(metrics$boundary_accuracy, 1)
  expect_equal(metrics$sparse_region_accuracy, 1)
})

test_that("adjusted Rand index handles perfect and degenerate partitions", {
  truth <- factor(c("A", "A", "B", "B"))
  crossed <- factor(c("A", "B", "A", "B"))
  one_cluster <- factor(rep("one", 4))

  perfect <- evaluate_spatial_refinement(truth, truth, truth)
  imperfect <- evaluate_spatial_refinement(truth, truth, crossed)
  degenerate <- evaluate_spatial_refinement(one_cluster, one_cluster, one_cluster)

  expect_equal(perfect$ari, 1)
  expect_equal(imperfect$ari, -0.5)
  expect_equal(degenerate$ari, 1)
  expect_true(is.na(perfect$correction_recall))
})

test_that("benchmark runner times named methods on identical inputs", {
  sim <- simulate_spatial_domains(
    n = 300, pattern = "jagged_stripes", noise = 0.20, seed = 401
  )
  methods <- list(
    identity = function(xy, labels) labels,
    oracle = function(xy, labels, samples) sim$truth
  )
  results <- benchmark_spatial_refiners(sim, methods, seed = 9)

  expect_identical(results$method, c("Initial", "identity", "oracle"))
  expect_equal(results$accuracy[results$method == "identity"],
               mean(sim$labels == sim$truth))
  expect_equal(results$accuracy[results$method == "oracle"], 1)
  expect_true(all(results$seconds >= 0))
  expect_true(all(is.na(results$error)))
})

test_that("benchmark runner can record method errors", {
  sim <- simulate_gradient_regions(n = 300, seed = 402)
  results <- benchmark_spatial_refiners(
    sim,
    list(fails = function(xy, labels) stop("deliberate failure")),
    include_initial = FALSE,
    on_error = "record"
  )

  expect_match(results$error, "deliberate failure")
  expect_true(is.na(results$accuracy))
  expect_equal(results$n, 300)
})

test_that("spatial benchmark constructor validates optional strata", {
  xy <- cbind(x = seq_len(5), y = seq_len(5))
  expect_error(
    spatial_benchmark(xy, letters[1:5], letters[1:5], boundary = 1:5),
    "logical"
  )
  expect_error(
    spatial_benchmark(xy, letters[1:4], letters[1:5]),
    "labels"
  )
})

test_that("licensed DLPFC and CRC benchmarks are bundled", {
  status <- available_spatial_benchmarks()
  expect_true(status$included[status$dataset == "dlpfc"])
  expect_false(status$included[status$dataset == "merfish"])
  expect_true(status$included[status$dataset == "crc"])

  dlpfc <- load_spatial_benchmark("dlpfc", scenario = 1)
  expect_s3_class(dlpfc, "spatial_refinement_benchmark")
  expect_equal(dim(dlpfc$xy), c(47329L, 2L))
  expect_equal(nlevels(dlpfc$truth), 7L)
  expect_equal(length(unique(dlpfc$samples)), 12L)
  expect_equal(mean(dlpfc$labels != dlpfc$truth), 2367 / 47329)
  expect_match(dlpfc$metadata$scenario_id, "DLPFC_")

  set.seed(17)
  old_seed <- .Random.seed
  crc <- load_spatial_benchmark(
    "crc", scenario = "CRC_random_25_r1"
  )
  expect_identical(.Random.seed, old_seed)
  expect_s3_class(crc, "spatial_refinement_benchmark")
  expect_equal(dim(crc$xy), c(194541L, 2L))
  expect_equal(nlevels(crc$truth), 19L)
  expect_equal(length(unique(crc$samples)), 1L)
  expect_equal(mean(crc$labels != crc$truth), 48635 / 194541)
  expect_silent(fibermargin:::.validate_class_anchors(
    crc$truth, crc$labels, crc$samples
  ))
  expect_identical(crc$metadata$license,
                   "Creative Commons Attribution 4.0 International (CC BY 4.0)")
  expect_identical(
    crc$labels,
    load_spatial_benchmark("crc", "CRC_random_25_r1")$labels
  )
  expect_identical(
    crc$labels,
    load_spatial_benchmark("colorectal", "CRC_random_25_r1")$labels
  )
  difficult_crc <- load_spatial_benchmark("crc", scenario = 59L)
  expect_silent(fibermargin:::.validate_class_anchors(
    difficult_crc$truth, difficult_crc$labels,
    difficult_crc$samples
  ))

  expect_error(load_spatial_benchmark("merfish"), "not bundled")
  expect_error(load_spatial_benchmark("dlpfc", 100), "Unknown")
  expect_error(load_spatial_benchmark("crc", "missing"), "Unknown")
})
