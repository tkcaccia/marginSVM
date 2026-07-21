test_that("categorical masks retain representation and diagnostics", {
  truth <- matrix(rep(c("background", "region"), each = 800), nrow = 40)
  truth[1:3, 1:3] <- NA_character_
  noisy <- corrupt_categorical_mask(truth, "impulse", rate = 0.12, seed = 2101)
  cleaned <- clean_categorical_mask(noisy, workers = 1)

  expect_identical(dim(cleaned), dim(truth))
  expect_identical(dimnames(cleaned), dimnames(truth))
  expect_type(cleaned, "character")
  expect_true(all(is.na(cleaned) == is.na(truth)))
  expect_identical(dim(attr(cleaned, "candidate")), dim(truth))
  expect_identical(dim(attr(cleaned, "repair_margin")), dim(truth))
  expect_identical(dim(attr(cleaned, "changed")), dim(truth))
  expect_true(all(is.na(attr(cleaned, "changed")) == is.na(truth)))
  expect_equal(
    attr(cleaned, "repair_margin"),
    attr(cleaned, "margin_score") - attr(cleaned, "required")
  )
  expect_null(attr(cleaned, "risk"))
  expect_null(attr(cleaned, "repair_score"))
})

test_that("factor and integer masks preserve their label representation", {
  factor_mask <- structure(
    factor(rep(c("A", "B"), each = 50), levels = c("unused", "A", "B")),
    dim = c(10L, 10L)
  )
  factor_cleaned <- clean_categorical_mask(factor_mask, workers = 1)
  expect_s3_class(factor_cleaned, "factor")
  expect_identical(levels(factor_cleaned), levels(factor_mask))
  expect_identical(dim(factor_cleaned), dim(factor_mask))

  integer_mask <- matrix(rep(c(10L, 20L), each = 50), nrow = 10)
  integer_cleaned <- clean_categorical_mask(integer_mask, workers = 1)
  expect_type(integer_cleaned, "integer")
  expect_true(all(integer_cleaned %in% c(10L, 20L)))
})

test_that("mask corruptions are deterministic and change the requested cells", {
  truth <- matrix(rep(c("A", "B", "C"), each = 400), nrow = 30)
  for (mechanism in c("impulse", "boundary", "patch")) {
    first <- corrupt_categorical_mask(truth, mechanism, rate = 0.20, seed = 2201)
    second <- corrupt_categorical_mask(truth, mechanism, rate = 0.20, seed = 2201)
    expect_identical(first, second)
    expect_equal(sum(first != truth), round(0.20 * length(truth)))
    expect_true(all(first %in% unique(as.vector(truth))))
    expect_silent(fibermargin:::.validate_class_anchors(
      as.vector(truth), as.vector(first)
    ))
  }
})

test_that("mask corruptions retain rare-class anchors", {
  truth <- matrix("A", nrow = 20L, ncol = 20L)
  truth[10L, 10L] <- "B"

  for (mechanism in c("impulse", "boundary", "patch")) {
    noisy <- corrupt_categorical_mask(truth, mechanism, rate = 0.50, seed = 2301)
    expect_equal(sum(noisy != truth), 200L)
    expect_silent(fibermargin:::.validate_class_anchors(
      as.vector(truth), as.vector(noisy)
    ))
  }
  expect_error(
    corrupt_categorical_mask(truth, "impulse", rate = 0.999, seed = 2302),
    "uncorrupted exemplar"
  )
})

test_that("mask metrics obey the correction-damage decomposition", {
  reference <- matrix(rep(c("A", "B"), each = 50), nrow = 10)
  initial <- reference
  initial[c(4, 14, 24, 34)] <- ifelse(initial[c(4, 14, 24, 34)] == "A", "B", "A")
  cleaned <- initial
  cleaned[c(4, 14, 50)] <- reference[c(4, 14, 50)]

  metrics <- evaluate_mask_cleaning(reference, initial, cleaned)
  expect_s3_class(metrics, "mask_cleaning_metrics")
  expect_equal(metrics$repair_identity_error, 0, tolerance = 1e-15)
  expect_equal(
    metrics$accuracy_gain,
    metrics$repaired_fraction - metrics$damaged_fraction,
    tolerance = 1e-15
  )
  expect_true(metrics$mean_iou >= 0 && metrics$mean_iou <= 1)
  expect_true(metrics$mean_boundary_iou >= 0 && metrics$mean_boundary_iou <= 1)
  expect_true(metrics$rare_class_iou >= 0 && metrics$rare_class_iou <= 1)

  perfect <- evaluate_mask_cleaning(reference, reference, reference)
  expect_equal(perfect$mean_iou, 1)
  expect_equal(perfect$macro_dice, 1)
  expect_equal(perfect$mean_boundary_iou, 1)
  expect_equal(perfect$boundary_accuracy, 1)
})
