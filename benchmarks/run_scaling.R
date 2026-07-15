#!/usr/bin/env Rscript

out_dir <- file.path("benchmarks", "results", "reviewer_response")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
quick <- identical(Sys.getenv("SPATIAL_REFINE_QUICK"), "1")
sizes <- if (quick) c(5000L, 60000L) else c(5000L, 20000L, 60000L, 150000L, 500000L)
replicates <- seq_len(if (quick) 1L else 3L)
design <- expand.grid(n = sizes, dimensions = c(2L, 3L), replicate = replicates)
results <- vector("list", nrow(design))

for (i in seq_len(nrow(design))) {
  row <- design[i, ]
  message("Scaling ", i, "/", nrow(design), ": n=", row$n, ", d=", row$dimensions)
  stdout <- tempfile()
  stderr <- tempfile()
  status <- system2(
    "/usr/bin/time",
    c("-l", file.path(R.home("bin"), "Rscript"), "benchmarks/scale_case.R",
      row$n, row$dimensions, 800000L + i),
    stdout = stdout,
    stderr = stderr
  )
  values <- strsplit(readLines(stdout, warn = FALSE)[1L], ",", fixed = TRUE)[[1L]]
  timing <- readLines(stderr, warn = FALSE)
  rss_line <- grep("maximum resident set size", timing, value = TRUE)
  rss <- if (length(rss_line)) as.numeric(sub("^\\s*([0-9]+).*$", "\\1", rss_line[1L])) else NA_real_
  results[[i]] <- data.frame(
    n = as.integer(values[1L]),
    dimensions = as.integer(values[2L]),
    seconds = as.numeric(values[3L]),
    accuracy = as.numeric(values[4L]),
    peak_rss_mb = rss / 1024^2,
    replicate = row$replicate,
    status = status
  )
  unlink(c(stdout, stderr))
}

metrics <- do.call(rbind, results)
write.csv(metrics, file.path(out_dir, "scaling_500k_metrics.csv"), row.names = FALSE)

suppressPackageStartupMessages(library(ggplot2))
p <- ggplot(metrics, aes(n, seconds, color = factor(dimensions))) +
  geom_point() +
  geom_line(aes(group = interaction(dimensions, replicate))) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Observations", y = "Elapsed seconds", color = "Dimensions") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "scaling_500k_runtime.png"), p, width = 7, height = 5, dpi = 220)

p_memory <- ggplot(metrics, aes(n, peak_rss_mb, color = factor(dimensions))) +
  geom_point() +
  geom_line(aes(group = interaction(dimensions, replicate))) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Observations", y = "Peak process memory (MB)", color = "Dimensions") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "scaling_500k_memory.png"), p_memory, width = 7, height = 5, dpi = 220)
