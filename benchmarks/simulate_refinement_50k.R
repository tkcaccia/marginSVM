#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(SpatialGraphRefine))

run_case <- function(dimensions, n = 60000L) {
  message("\n--- ", dimensions, "D benchmark with n = ", n, " ---")
  sim <- simulate_spatial_clusters(
    n = n,
    dimensions = dimensions,
    k = 6,
    samples = 2,
    noise = 0.12,
    seed = 100 + dimensions
  )

  tiles <- rep(if (dimensions == 2L) 8L else 5L, dimensions)
  elapsed <- system.time({
    refined <- refine_spatial_svm(
      sim$xy,
      sim$labels,
      samples = sim$samples,
      tiles = tiles,
      n_features = 96,
      epochs = 6,
      gamma = 0.35,
      seed = 42
    )
  })

  before <- mean(sim$labels == sim$truth)
  after <- mean(refined == sim$truth)
  print(elapsed)
  cat("accuracy_before:", signif(before, 4), "\n")
  cat("accuracy_after: ", signif(after, 4), "\n")
  invisible(refined)
}

run_case(2L)
run_case(3L)
