#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(marginSVM))
args <- commandArgs(trailingOnly = TRUE)
n <- as.integer(args[1L])
dimensions <- as.integer(args[2L])
workers <- as.integer(args[3L])
seed <- as.integer(args[4L])
pattern <- if (dimensions == 3L) "layers3d" else "jagged_stripes"
sim <- simulate_spatial_domains(
  n, pattern, dimensions = dimensions, noise = 0.25,
  density_profile = "uniform", seed = seed
)
invisible(gc())
timing <- system.time(refined <- marginSVM:::.refine_spatial_svm_engine(
  sim$xy, sim$labels, sim$samples,
  control = list(workers = workers, seed = seed + 100000L)
))
cat(paste(n, dimensions, workers, unname(timing["elapsed"]),
          mean(refined == sim$truth), attr(refined, "tiles"), sep = ","), "\n")
