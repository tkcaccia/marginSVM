#!/usr/bin/env Rscript

out_dir <- file.path("benchmarks", "results", "marginsvm_scaling")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sizes <- c(5000L, 20000L, 60000L, 150000L, 500000L)
design <- expand.grid(n = sizes, dimensions = c(2L, 3L), replicate = 1:3)
design <- design[design$n < 500000L | design$replicate == 1L, ]
results <- vector("list", nrow(design))

for (i in seq_len(nrow(design))) {
  d <- design[i, ]
  message("marginSVM scaling ", i, "/", nrow(design), ": n=", d$n,
          ", d=", d$dimensions)
  stdout <- tempfile()
  stderr <- tempfile()
  status <- system2(
    "/usr/bin/time",
    c("-l", file.path(R.home("bin"), "Rscript"),
      "benchmarks/marginsvm_scale_case.R", d$n, d$dimensions, 4L, 2500000L + i),
    stdout = stdout, stderr = stderr
  )
  values <- strsplit(readLines(stdout, warn = FALSE)[1L], ",", fixed = TRUE)[[1L]]
  timing <- readLines(stderr, warn = FALSE)
  rss_line <- grep("maximum resident set size", timing, value = TRUE)
  rss <- if (length(rss_line)) as.numeric(sub("^\\s*([0-9]+).*$", "\\1", rss_line[1L])) else NA_real_
  results[[i]] <- data.frame(
    n = as.integer(values[1L]), dimensions = as.integer(values[2L]),
    workers = as.integer(values[3L]), seconds = as.numeric(values[4L]),
    accuracy = as.numeric(values[5L]), tiles = as.integer(values[6L]),
    peak_rss_mb = rss / 1024^2, replicate = d$replicate, status = status
  )
  unlink(c(stdout, stderr))
}

metrics <- do.call(rbind, results)
write.csv(metrics, file.path(out_dir, "marginsvm_scaling_metrics.csv"), row.names = FALSE)

suppressPackageStartupMessages(library(ggplot2))
runtime <- ggplot(metrics, aes(n, seconds, color = factor(dimensions))) +
  geom_point() + geom_line(aes(group = interaction(dimensions, replicate))) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Observations", y = "Elapsed seconds", color = "Dimensions") +
  theme_bw(base_size = 10)
ggsave(file.path(out_dir, "marginsvm_scaling_runtime.png"), runtime,
       width = 7.2, height = 4.8, dpi = 240, bg = "white")

memory <- ggplot(metrics, aes(n, peak_rss_mb, color = factor(dimensions))) +
  geom_point() + geom_line(aes(group = interaction(dimensions, replicate))) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Observations", y = "Peak process memory (MB)", color = "Dimensions") +
  theme_bw(base_size = 10)
ggsave(file.path(out_dir, "marginsvm_scaling_memory.png"), memory,
       width = 7.2, height = 4.8, dpi = 240, bg = "white")

print(aggregate(cbind(seconds, peak_rss_mb, accuracy) ~ n + dimensions,
                metrics, mean))
