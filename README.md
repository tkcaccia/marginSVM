# marginSVM

[![R CMD check](https://github.com/tkcaccia/marginSVM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tkcaccia/marginSVM/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/tkcaccia/marginSVM/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/tkcaccia/marginSVM/actions/workflows/pkgdown.yaml)

`marginSVM` refines noisy spatial cluster assignments without rerunning the
upstream clustering method. It uses coordinates and labels, with an optional
tissue identifier for datasets containing multiple sections. The C++ engine
combines overlapping local Nystrom SVMs, graph evidence, and edge-aware
total-variation decoding.

## Installation

```r
install.packages("remotes")
remotes::install_github("tkcaccia/marginSVM", build_vignettes = TRUE)
```

## Basic use

```r
library(marginSVM)

sim <- simulate_gradient_regions(
  n = 50000,
  minority = 0.25,
  samples = 2,
  density_profile = "strong"
)

refined <- refine_spatial_svm(
  sim$xy,
  sim$labels,
  samples = sim$samples
)

mean(sim$labels == sim$truth)
mean(refined == sim$truth)
```

Only `xy` and `labels` are required. Tiling, overlap, local models, graph fusion,
and decoding are automatic. Each value of `samples` is processed independently.

## Large datasets and backends

Large 2D and 3D inputs are automatically divided into overlapping tiles. CPU
tile models use up to four physical cores by default; the worker count can be
limited directly:

```r
refined <- refine_spatial_svm(xy, labels, samples, workers = 8)
```

The built-in backend is C++ CPU. Optional CUDA or Metal providers can implement
the same score contract:

```r
spatial_svm_backend_capabilities()
refined <- refine_spatial_svm(xy, labels, samples, backend = "auto")
```

No native GPU performance claim is made by this release.

## Diagnostics

The returned factor preserves the original label levels. Compact diagnostics are
stored as attributes:

```r
attr(refined, "confidence")
attr(refined, "margin")
attr(refined, "local_support")
attr(refined, "tiles")
attr(refined, "abstained_samples")
```

See the package vignette for visual examples and a Seurat workflow:

```r
vignette("marginSVM", package = "marginSVM")
```

## Validation

The manuscript evaluates 708 simulations spanning eleven geometries, 2D/3D
layered gradients, five region-density profiles, four error mechanisms, and one
or three tissues. It also reports controlled label-corruption experiments on
194,541 measured locations and 19 biological annotation classes from a real
VisiumHD colorectal tissue.

The complete benchmark and manuscript statistics are regenerated with:

```sh
Rscript benchmarks/marginsvm_complete_simulation_benchmark.R
Rscript benchmarks/analyze_complete_simulation_for_manuscript.R
```

The [project site](https://tkcaccia.github.io/marginSVM/) contains the function
reference, vignettes, benchmark summary, and manuscript.
