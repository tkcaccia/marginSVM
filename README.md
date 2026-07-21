# FiberMargin

[![R CMD check](https://github.com/tkcaccia/fibermargin/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tkcaccia/fibermargin/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/tkcaccia/fibermargin/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/tkcaccia/fibermargin/actions/workflows/pkgdown.yaml)

FiberMargin performs training-free repair of categorical masks and spatial labels.
It is a post-processing method: upstream segmentation, clustering, or annotation
remains unchanged. A deterministic atlas of rotated Hilbert paths transports
leave-self-out class evidence from both directions using normalized path length
only. Its class score is the two-sided enclosure
`(sqrt(left) + sqrt(right))^2`. FiberMargin changes a label only when its
full-atlas margin clears a conservative chart-disagreement barrier; locally
isolated points receive extra protection at that decision.

The package exposes one production refiner with four arguments. It does not need
expression values, a graph, training labels, or user-selected tuning parameters.

## Installation

```r
install.packages("remotes")
remotes::install_github("tkcaccia/fibermargin", build_vignettes = TRUE)
```

## Basic Use

```r
library(fibermargin)

sim <- simulate_gradient_regions(
  n = 50000,
  minority = 0.25,
  samples = 2,
  density_profile = "strong"
)

refined <- refine_spatial_labels(
  xy = sim$xy,
  labels = sim$labels,
  samples = sim$samples
)

c(
  input = mean(sim$labels == sim$truth),
  refined = mean(refined == sim$truth)
)
```

Only `xy` and `labels` are required. When `samples` is supplied, every specimen
is refined independently. `workers` is a CPU budget shared across independent
specimens and independent atlas charts:

```r
refined <- refine_spatial_labels(xy, labels, samples, workers = 4)
```

Chart parallelism is deterministic: changing `workers` does not change labels or
diagnostics. In the package benchmark, four workers provided 2.01--2.52x speedups
for single 100,000-point 2D/3D fields and 3.37x across four specimens.

## Clean a Categorical Mask

For a 2D pixel mask or 3D voxel array, coordinates are generated automatically.
Missing cells remain missing. No source image, intensity, model logit, or clean-mask
training collection is used.

```r
truth <- matrix(rep(c("background", "region"), each = 5000), nrow = 100)
noisy <- corrupt_categorical_mask(
  truth, mechanism = "boundary", rate = 0.15, seed = 4
)
cleaned <- clean_categorical_mask(noisy, workers = 1)

evaluate_mask_cleaning(truth, noisy, cleaned)
```

The evaluator reports multiclass mean IoU, one-grid-step Boundary IoU, rare-class
IoU, correction recall, and damage. Its accuracy obeys the exact identity
`gain = repaired_fraction - damaged_fraction`.

## Diagnostics

The result is an ordinary factor with the input levels. Pointwise diagnostics are
attached without expanding the public interface:

```r
attr(refined, "candidate")
attr(refined, "margin_score")
attr(refined, "required")
attr(refined, "repair_margin")
attr(refined, "atlas_dispersion")
attr(refined, "isolation")
attr(refined, "changed")
```

`margin_score` is an uncalibrated rival-versus-observed evidence contrast, not a
probability or expected loss. `repair_margin >= 0` is the selective change
criterion for a candidate that differs from the observed label.
`atlas_dispersion` is a deterministic chart-disagreement scale, not a standard
error or calibrated probability.
`isolation` is the multiclass local-gap protection factor; it is one for the
binary ballot.

## Evaluation

The package provides planar, gradient, complex-shape, and volumetric simulators,
plus one evaluator used by every benchmark:

```r
evaluate_spatial_refinement(
  truth = sim$truth,
  initial = sim$labels,
  refined = refined,
  boundary = sim$boundary,
  regions = sim$area,
  sparse = sim$sparse,
  method = "FiberMargin"
)
```

`correction_recall` is the fraction of initially wrong labels repaired.
`damage_rate` is the fraction of initially correct labels made wrong. Worst-class,
sparse-region, and boundary accuracy prevent aggregate accuracy from hiding local
failures.

Compact DLPFC and CRC Visium HD coordinate/annotation benchmarks are bundled. CRC
denotes colorectal cancer. DLPFC is redistributed under the `spatialLIBD` terms. The
CRC derivative is CC BY 4.0 and contains only coordinates, 19 author-derived region labels, and stored
corruption recipes; expression counts and tissue imagery are excluded. MERFISH
inputs remain external. See
`inst/extdata/REAL_DATA_LICENSES.md` for sources, attribution, and change notices.

```r
crc <- load_spatial_benchmark("crc", "CRC_random_25_r1")
mean(crc$labels != crc$truth)
```

## Current Validation

The C++ implementation is evaluated on predefined CRC and MERFISH corruption
protocols. On the direct-comparator matrices, it obtains 0.8600 accuracy and
0.7812 ARI across 60 CRC corruptions, within 0.0005 and 0.0016 of multiscale
mode, respectively. Across 45 MERFISH corruptions it obtains 0.8886 accuracy
and 0.7951 ARI, exceeding the strongest fixed comparator by 0.0044 accuracy
and 0.0063 ARI. These are controlled label-recovery results, not claims that
the biological annotations are naturally wrong.

The protocol definitions, package-build metadata, per-case ledgers, and matched plots
are under `benchmarks/results/fibermargin_final_enclosure_real_matrix_v2/`.
Additional simulation and mask panels are retained in the repository and are
re-evaluated only when the production source changes.

See `vignette("fibermargin", package = "fibermargin")` for the workflow,
`vignette("benchmarks", package = "fibermargin")` for benchmark design,
and `vignette("fibermargin_reference", package = "fibermargin")` for the package
reference API.
