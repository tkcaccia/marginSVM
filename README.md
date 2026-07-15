# marginSVM

[![R CMD check](https://github.com/tkcaccia/marginSVM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tkcaccia/marginSVM/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/tkcaccia/marginSVM/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/tkcaccia/marginSVM/actions/workflows/pkgdown.yaml)

The `SpatialGraphRefine` R package provides two complementary post-processors for noisy 2D and
3D tissue-domain assignments. `refine_spatial_svm()` implements marginSVM
(Multiscale Adaptive Refinement with Graph-Integrated Nyström SVM), built around
a structured C++ multiclass engine with cross-fitted Nyström
features, overlapping tiles, multiple tissues, and multicore execution.
`refine_spatial_clusters()` retains
the conservative kNN-graph alternative for locally isolated errors.

The usual workflow has only two required inputs:

## Install locally

```r
install.packages("Rcpp")
remotes::install_github("tkcaccia/marginSVM")
```

For development from a local checkout:

```sh
R CMD INSTALL .
```

## Example

```r
library(SpatialGraphRefine)

sim <- simulate_gradient_regions(
  50000, minority = 0.25, samples = 2, density_profile = "strong"
)
refined <- refine_spatial_svm(sim$xy, sim$labels, sim$samples)

mean(sim$labels == sim$truth)
mean(refined == sim$truth)
table(sim$area, refined)
```

## Large tiled datasets

Large SVM inputs are partitioned automatically into prediction cores with
overlapping training halos:

```r
refined <- refine_spatial_svm(xy, labels, samples)
```

Every tissue is partitioned independently. Automatic execution combines regular
spatial tiles with adaptive boundary-rich tiles targeting about 5,000 core
observations. Cross-fitted local models return continuous multiclass scores, which
are blended across overlapping halos and reconciled by edge-aware total variation.
High-cardinality annotations also use topology-aware abstention for coherent label
regions that cannot be identified as errors from coordinates alone.

An opt-in high-class extension adds anisotropic multiscale neighborhoods,
pairwise boundary SVMs, a continuous trust field, pointwise
retain/change/unresolved decisions, and stability-based rare-region protection:

```r
refined_v2 <- refine_spatial_svm(
  xy, labels, samples,
  control = list(experimental_v2 = 1)
)
attr(refined_v2, "decision")
attr(refined_v2, "trust")
```

The extension improved the frozen 19-class colorectal confirmation but not a
mixed 2D/3D geometry block, so v1 remains the general default.

```r
control <- list(
  overlap = 0.25,
  target_tile_size = 5000,
  workers = 6
)
refined <- refine_spatial_svm(xy, labels, samples, control = control)
```

The default 25% halo, 48-landmark, 8-epoch, 24-iteration profile was frozen after
paired tuning and checked on 560 independent simulated scenarios. It was 1.68-fold
faster overall and 1.80-fold faster in the 19-class subset than the previous
40%/64/10/40 profile. Four CPU workers processed 500,000 points in 6.19 seconds in
2D and 10.14 seconds in 3D on the benchmark machine.

`spatial_svm_backend_capabilities()` reports available structured-SVM backends.
CPU is built in. CUDA and Metal implementations can register through
`register_spatial_svm_backend()`. An unavailable requested accelerator warns and
uses CPU. Tile fitting uses a portable native C++ thread pool.

## Benchmark

The complete simulation matrix covers all eleven geometric structures, layered
A-B-B-C gradients, five region-concentration profiles, four error mechanisms,
2D/3D coordinates, and single/multiple tissues:

```sh
Rscript benchmarks/marginsvm_complete_simulation_benchmark.R
Rscript benchmarks/analyze_complete_simulation_for_manuscript.R
```

For the frozen synthetic study and manuscript figures:

```sh
Rscript benchmarks/reviewer_response_experiments.R
Rscript benchmarks/svm_gradient_benchmark.R
Rscript benchmarks/reviewer_additional_experiments.R
Rscript benchmarks/run_scaling.R
Rscript benchmarks/summarize_manuscript_results.R
```

These scripts write raw CSV metrics and publication PNG figures to
`benchmarks/results/reviewer_response/`, including:

- 1,120 held-out combinations of ten geometries, four noise mechanisms, four
  corruption fractions, and seven confirmatory seeds;
- published SpaGCN and GraphST label refiners, C++ kNN, Potts-like ICM, and
  image-morphology comparators;
- rare-domain, thin-layer, boundary, anisotropy, duplicate-coordinate, and
  irregular-density stress tests;
- rotation, translation, row-order, label-permutation, and multi-sample checks;
- CPU runtime and peak memory through 500,000 observations in 2D and 3D.
- layered A-B-B-C gradient mixtures at 5--40%, in 2D/3D and one/three tissues.
- a held-out multiscale SVM study across gradient mixtures and five nonlinear
  geometric families, including visual border comparisons.

The manuscript includes disjoint synthetic seeds and a third untouched
60-scenario seed block for the 19-class colorectal VisiumHD stress test.
Independent biological datasets remain reserved for external validation.

See `docs/literature_benchmark_notes.md` for related tools and benchmark papers.

## Implementation scope

The graph refiner uses an exact CPU kd-tree. The structured SVM uses boundary-
enriched Nyström features, robust two-fold multiclass hinge optimization,
uncertainty-gated graph evidence, and primal-dual TV decoding in C++. CUDA and
Metal providers are supported by contract; this release makes no GPU performance
claim.
