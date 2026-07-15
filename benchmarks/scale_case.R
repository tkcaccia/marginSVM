#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(SpatialGraphRefine))
args <- commandArgs(trailingOnly = TRUE)
n <- as.integer(args[1L])
dimensions <- as.integer(args[2L])
seed <- as.integer(args[3L])
pattern <- if (dimensions == 3L) "layers3d" else "jagged_stripes"
sim <- simulate_spatial_domains(n, pattern, dimensions = dimensions, noise = 0.25, seed = seed)
invisible(gc())
timing <- system.time(refined <- refine_spatial_clusters(sim$xy, sim$labels))
cat(paste(n, dimensions, unname(timing["elapsed"]), mean(refined == sim$truth), sep = ","), "\n")
