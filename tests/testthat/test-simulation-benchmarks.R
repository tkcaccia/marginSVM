test_that("all simulators return the common benchmark fields", {
  simulations <- list(
    clusters = simulate_spatial_clusters(n = 300, noise = 0.20, seed = 501),
    gradient = simulate_gradient_regions(n = 300, minority = 0.05, seed = 502),
    domains = simulate_spatial_domains(
      n = 300, pattern = "rings", noise = 0.20, seed = 503
    ),
    complex = simulate_complex_spatial_domains(
      n = 300, shape = "tubular_network", noise = 0.20, seed = 504
    ),
    volume = simulate_volumetric_domains(
      n = 300, shape = "folded_layers", noise = 0.20, seed = 505
    )
  )

  for (simulation in simulations) {
    expect_s3_class(simulation, "spatial_refinement_benchmark")
    expect_equal(nrow(simulation$xy), 300L)
    expect_length(simulation$labels, 300L)
    expect_length(simulation$truth, 300L)
    expect_length(simulation$samples, 300L)
    expect_type(simulation$boundary, "logical")
    expect_type(simulation$sparse, "logical")
    expect_equal(sum(simulation$boundary), 60L, tolerance = 2L)
    expect_equal(simulation$corrupted, simulation$labels != simulation$truth)
    expect_silent(fibermargin:::.validate_class_anchors(
      simulation$truth, simulation$labels, simulation$samples
    ))
  }
  expect_equal(ncol(simulations$volume$xy), 3L)
})

test_that("complex and volumetric geometry catalogues are reproducible", {
  complex_a <- simulate_complex_spatial_domains(
    n = 300, shape = "3D shells", density_profile = "extreme",
    noise_type = "patch", samples = 2, seed = 510
  )
  complex_b <- simulate_complex_spatial_domains(
    n = 300, shape = "shells_3d", density_profile = "extreme",
    noise_type = "patch", samples = 2, seed = 510
  )
  volume_a <- simulate_volumetric_domains(
    n = 300, shape = "Helical channels", acquisition = "irregular_z",
    noise_type = "region", samples = 2, seed = 511
  )
  volume_b <- simulate_volumetric_domains(
    n = 300, shape = "helical_channels", acquisition = "irregular_z",
    noise_type = "region", samples = 2, seed = 511
  )

  expect_identical(complex_a$xy, complex_b$xy)
  expect_identical(complex_a$labels, complex_b$labels)
  expect_identical(volume_a$xy, volume_b$xy)
  expect_identical(volume_a$labels, volume_b$labels)
  expect_equal(length(unique(volume_a$samples)), 2L)
})

test_that("cluster simulator corrupts the requested number without no-op labels", {
  sim <- simulate_spatial_clusters(n = 500, k = 4, noise = 0.20, seed = 520)
  expect_equal(sum(sim$corrupted), 100L)
  expect_equal(mean(sim$labels != sim$truth), 0.20)
  high_noise <- simulate_spatial_clusters(n = 500, k = 4, noise = 0.98, seed = 521)
  expect_equal(sum(high_noise$corrupted), 490L)
  expect_silent(fibermargin:::.validate_class_anchors(
    high_noise$truth, high_noise$labels, high_noise$samples
  ))
  expect_error(simulate_spatial_clusters(n = 100, noise = 1), "noise")
})
