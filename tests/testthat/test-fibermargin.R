test_that("FiberMargin returns stable factors and diagnostics", {
  simulation <- simulate_spatial_domains(
    n = 1200, pattern = "jagged_stripes", noise = 0.20,
    noise_type = "boundary", seed = 1501
  )
  refined <- refine_spatial_labels(
    simulation$xy, simulation$labels, simulation$samples, workers = 1
  )

  expect_s3_class(refined, "factor")
  expect_identical(levels(refined), levels(factor(simulation$labels)))
  expect_length(refined, nrow(simulation$xy))
  expect_length(attr(refined, "candidate"), nrow(simulation$xy))
  expect_true(all(is.finite(attr(refined, "margin_score"))))
  expect_true(all(is.finite(attr(refined, "required"))))
  expect_equal(
    attr(refined, "repair_margin"),
    attr(refined, "margin_score") - attr(refined, "required")
  )
  expect_null(attr(refined, "risk"))
  expect_true(all(attr(refined, "atlas_dispersion") >= 0))
  expect_identical(attr(refined, "changed"), refined != factor(simulation$labels))

  observed <- as.integer(factor(simulation$labels, levels = levels(refined)))
  candidate <- as.integer(attr(refined, "candidate"))
  expect_true(all(attr(refined, "required") >= 0))
  expect_identical(
    attr(refined, "changed"),
    candidate != observed & attr(refined, "repair_margin") >= 0
  )
})

test_that("binary local-volume evidence uses the shared margin admission rule", {
  set.seed(1551)
  xy <- matrix(runif(1600), ncol = 2)
  labels <- factor(ifelse(xy[, 1] + 0.12 * sin(8 * xy[, 2]) < 0.5, "A", "B"))
  noisy <- labels
  corrupted <- sample.int(length(noisy), 120L)
  noisy[corrupted] <- ifelse(noisy[corrupted] == "A", "B", "A")
  noisy <- factor(noisy, levels = levels(labels))

  refined <- refine_spatial_labels(xy, noisy, workers = 1)
  observed <- as.integer(noisy)
  candidate <- as.integer(attr(refined, "candidate"))

  expect_true(all(attr(refined, "required") == 1))
  expect_true(all(attr(refined, "margin_score")[candidate == observed] == 0))
  expect_true(all(candidate[attr(refined, "margin_score") == 0] ==
                  observed[attr(refined, "margin_score") == 0]))
  expect_identical(
    attr(refined, "changed"),
    candidate != observed & attr(refined, "repair_margin") >= 0
  )
})

test_that("specimens are refined independently", {
  first <- simulate_spatial_domains(
    n = 800, pattern = "rings", noise = 0.25, seed = 1601
  )
  second <- simulate_spatial_domains(
    n = 900, pattern = "branching", noise = 0.25, seed = 1602
  )
  xy <- rbind(first$xy, second$xy + 20)
  labels <- factor(
    c(as.character(first$labels), as.character(second$labels)),
    levels = union(levels(first$labels), levels(second$labels))
  )
  samples <- factor(c(rep("first", 800), rep("second", 900)))

  joint <- refine_spatial_labels(xy, labels, samples, workers = 1)
  separate_first <- refine_spatial_labels(
    first$xy, factor(first$labels, levels = levels(labels)), workers = 1
  )
  separate_second <- refine_spatial_labels(
    second$xy, factor(second$labels, levels = levels(labels)), workers = 1
  )

  expect_identical(as.character(joint[seq_len(800)]), as.character(separate_first))
  expect_identical(as.character(joint[800 + seq_len(900)]), as.character(separate_second))
})

test_that("chart and specimen workers are bitwise deterministic", {
  simulation <- simulate_spatial_domains(
    n = 2400, pattern = "jagged_stripes", noise = 0.25,
    samples = 3, seed = 1651
  )
  serial <- refine_spatial_labels(
    simulation$xy, simulation$labels, simulation$samples, workers = 1
  )
  parallel <- refine_spatial_labels(
    simulation$xy, simulation$labels, simulation$samples, workers = 2
  )

  expect_identical(attr(serial, "workers"), 1L)
  expect_identical(attr(parallel, "workers"), 2L)
  attr(serial, "workers") <- NULL
  attr(parallel, "workers") <- NULL
  expect_identical(parallel, serial)

  single <- simulate_spatial_domains(
    n = 1600, pattern = "rings", noise = 0.25, samples = 1, seed = 1652
  )
  serial_chart <- refine_spatial_labels(
    single$xy, single$labels, single$samples, workers = 1
  )
  parallel_chart <- refine_spatial_labels(
    single$xy, single$labels, single$samples, workers = 2
  )
  attr(serial_chart, "workers") <- NULL
  attr(parallel_chart, "workers") <- NULL
  expect_identical(parallel_chart, serial_chart)
})

test_that("zero-evidence multiclass ties retain the observed label", {
  xy <- matrix(0, nrow = 30, ncol = 2)
  labels <- factor(rep(c("A", "B", "C"), length.out = 30))

  refined <- refine_spatial_labels(xy, labels, workers = 1)

  expect_identical(as.character(refined), as.character(labels))
  expect_identical(
    as.character(attr(refined, "candidate")), as.character(labels)
  )
  expect_true(all(attr(refined, "margin_score") == 0))
  expect_true(all(attr(refined, "required") == 0))
})

test_that("an unused factor level is not introduced in a regular spatial sample", {
  set.seed(1701)
  xy <- matrix(runif(800), ncol = 2)
  labels <- factor(ifelse(xy[, 1] < 0.5, "A", "B"), levels = c("C", "A", "B"))
  refined <- refine_spatial_labels(xy, labels, workers = 1)
  expect_false(any(refined == "C"))
  expect_false(any(attr(refined, "candidate") == "C"))
})

test_that("FiberMargin exposes no mathematical tuning controls", {
  set.seed(1801)
  xy <- matrix(runif(82), ncol = 2)
  labels <- factor(c(rep("A", 20), rep("B", 21)), levels = c("A", "B"))
  expect_error(
    fibermargin:::.fiber_margin_engine(
      xy, labels, workers = 1,
      control = list(tail_fraction = 0.50)
    ),
    "no tuning controls"
  )
})

test_that("fixed-k modal reference uses the shared neighbour-mode kernel", {
  xy <- as.matrix(expand.grid(x = seq_len(5L), y = seq_len(5L)))
  labels <- factor(ifelse(xy[, 1L] <= 3L, "A", "B"))
  samples <- factor(rep("mask", nrow(xy)))

  modal <- fibermargin:::.local_modal_filter_labels(
    xy, labels, samples, neighbors = 8L
  )
  direct <- fibermargin:::.refine_published_labels(
    xy, labels, samples, method = "graphst", neighbors = 8L
  )

  expect_identical(modal, direct)
  expect_identical(attr(modal, "neighbors"), attr(direct, "neighbors"))
})

test_that("alpha-expansion Potts control preserves factor labels and diagnostics", {
  simulation <- simulate_spatial_domains(
    n = 800, pattern = "rings", noise = 0.20,
    noise_type = "random", seed = 1851
  )
  refined <- fibermargin:::.alpha_expansion_potts_labels(
    simulation$xy, simulation$labels, simulation$samples
  )

  expect_s3_class(refined, "factor")
  expect_identical(levels(refined), levels(simulation$labels))
  expect_length(refined, nrow(simulation$xy))
  expect_true(all(!is.na(refined)))
  expect_true(is.finite(attr(refined, "energy")))
  expect_true(attr(refined, "energy") >= 0)
  expect_true(all(attr(refined, "neighbors") >= 0L))
  expect_true(all(attr(refined, "changes") >= 0L))
})
