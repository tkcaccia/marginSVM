# Publication experiment plan

## Central claim

The method is a fast, clustering-agnostic post-processor that repairs spatially
inconsistent tissue-domain labels while retaining thin, curved, rare, and
disconnected structures. It is not an ab initio transcriptomic clustering method.

## Primary hypotheses

1. Refinement improves agreement with known domains over unrefined labels across
   spatial platforms and upstream clustering methods.
2. Discordance-gated graph refinement provides a better correction-preservation
   tradeoff than published direct label-refinement rules and local smoothing.
3. Runtime and memory scale to at least 500,000 observations on CPU.
4. Frozen defaults generalize across slides without dataset-specific tuning.

## Synthetic experiments: completed

- Geometries: jagged stripes, wavy and thin layers, rings, spirals, branching
  sectors, lobes, rare islands, disconnected domains, intermixed domains, and
  curved 3D layers.
- Sample sizes: 5k, 20k, 60k, 150k, and 500k observations.
- Domain counts: 2, 3, 5, 8, and 12.
- Noise fractions: 0%, 5%, 15%, 25%, and 40% where applicable.
- Noise mechanisms: independent flips, boundary errors, correlated patches, and
  coherent regional relabeling.
- Replication: three development seeds and seven untouched confirmatory seeds
  in the full factorial experiment.

The confirmatory experiment contains 1,120 datasets and 13.44 million unique
simulated observations. Endpoints include ARI, overall accuracy, minimum-domain
recall, analytic-boundary accuracy, changed-label precision, correction recall,
erroneous-change rate, runtime, and peak memory.

Additional completed studies cover clean-label damage, rare and thin structures,
domain count, ablations, parameter sensitivity, 3D coordinate anisotropy,
duplicated coordinates, irregular density and holes, regular-grid morphology,
invariance, and multi-sample stratification. CPU scaling reaches 500,000
observations in 2D and 3D.

## Real-data experiments: awaiting datasets

Use datasets with independent annotations and multiple technologies:

- 12 human DLPFC 10x Visium sections with cortical-layer annotations.
- Mouse hypothalamic MERFISH sections with annotated anatomical regions.
- Mouse brain Stereo-seq or Slide-seq data for high-resolution scalability.
- Breast or colorectal tumor Visium data with pathologist-defined regions.
- A multi-section or reconstructed 3D dataset for sample-stratified refinement.

For each dataset, generate several sets of initial labels on fixed representations.
These upstream procedures create test inputs and are not performance comparators.
Apply every direct post-processor to each identical label set to isolate refinement
quality from upstream feature learning.

## Comparators: synthetic phase completed

- No post-processing.
- SpaGCN's published strict-majority label-refinement function.
- GraphST's published 50-neighbor modal label-refinement function.
- Spatial k-nearest-neighbor majority vote.
- Potts-like fixed-graph iterative smoothing.
- Morphological modal filtering on regular grids.

End-to-end spatial clustering pipelines and the precursor KODAMA SVM are excluded
from method comparisons. For real data, run documented refiner defaults first,
tune only on validation slides, and report default-versus-tuned results separately.

## Statistical analysis

- Pre-specify ARI as the primary real-data endpoint and macro recall as the
  primary synthetic endpoint for any additional simulations.
- Compare paired outputs using a hierarchical bootstrap over subjects, slides,
  upstream clusterers, and seeds; report paired differences and 95% intervals.
- Control false discovery rate for secondary pairwise comparisons.
- Report erased rare regions, over-smoothed borders, and coherent-error failures.
- Do not treat spots within a slide as independent biological replicates.

## Remaining publication-critical work

1. Freeze the package and register the real-data protocol before opening the
   supplied datasets.
2. Run multiple technologies, tissue types, sections, and upstream clusterers;
   keep all post-processing paired to identical initial labels.
3. Add expert or pathology review of changed labels and boundary-specific errors.
4. Compare defaults and validation-only tuning with a hierarchical bootstrap at
   the subject or slide level.
5. Report abstention behavior, support calibration, and every observed failure;
   coherent regional errors are not identifiable from spatial labels alone.

## Reproducibility and release

Freeze package version, compiler, operating system, CPU, dependency versions,
random seeds, and dataset checksums. Publish scripts, raw metrics, session info,
and figure-generation code. Keep the real-data test set held out and do not
change defaults after inspecting it.
