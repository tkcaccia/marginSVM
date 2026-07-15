# Synthetic reviewer-response ledger

This document records how the simulation-only revision addresses the major
methodological concerns. Independently annotated biological data remain the one
publication-critical evidence class not yet available.

| Concern | Experiment or change | Outcome |
|:--|:--|:--|
| Method needed a scalable implementation | Implemented exact per-sample C++ kd-tree graph refinement | 500k observations required 3.18 s and 434 MB in 2D, or 4.77 s and 472 MB in 3D |
| Large datasets may exceed one-process memory | Added deterministic 2D/3D tile cores, overlapping halos, automatic tile sizing, and forked CPU workers | At 150k observations, 15% halos reproduced untiled labels exactly in 2D and at 99.9993% in 3D; four workers were 2.24x and 2.02x faster than one tiled worker |
| GPU requirements vary by workstation | Added backend capability reporting, strict/fallback selection, and a validated provider registration contract for CUDA and Metal | CPU is built in and tested; this release makes no unimplemented GPU performance claim |
| Excess user tuning | Reduced the main interface to coordinates, labels, and optional sample identifiers; froze an automatic sample-specific neighborhood rule | All six simulated samples independently selected 20 neighbors |
| Risk of boundary and rare-domain erosion | Added current-label discordance, consensus, margin, and preservation gates; tested analytic boundaries, rare islands, and thin layers | Boundary accuracy exceeded matched kNN by 0.0308 and clean-label damage was lower by 0.0186; a thin-layer resolution limit was quantified |
| Benchmark was narrow | Added ten 2D geometries, curved 3D layers, 2--12 domains, irregular density, holes, duplicated coordinates, anisotropy, and four error mechanisms | Confirmatory matrix contains 1,120 held-out datasets and 13.44 million unique observations |
| No strong direct comparators | Added source-faithful SpaGCN and GraphST label refiners, matched compiled kNN, Potts-like ICM, and regular-grid morphology; excluded end-to-end clustering and the precursor KODAMA SVM | GraphST was more accurate but damaged correct labels much more often; SpatialGraphRefine outperformed SpaGCN on accuracy, boundaries, and preservation |
| Development leakage | Used seeds 1--3 for development, froze defaults, and reported only seeds 4--10 in the primary factorial study | Paired bootstrap estimates use matched held-out seeds only |
| Failure modes were hidden | Added coherent region swaps, correlated patches, clean labels, extreme duplicates, sub-neighborhood layers, and uncorrected 3D units | Coherent swaps were not recoverable; patch errors changed little; extreme duplicates, thin layers, and mismatched coordinate units degraded performance |
| Diagnostic score could be mistaken for probability | Renamed the diagnostic to local support and measured calibration | ECE was 0.115; the manuscript explicitly rejects a posterior-probability interpretation |
| Reproducibility was incomplete | Added raw metrics, frozen scripts, unit tests, source package build, machine/compiler details, and figure-generation scripts | 90 test expectations pass and `R CMD check --no-manual` reports `Status: OK` |
| Biological validity was overstated | Reframed the contribution as a conservative label corrector and separated synthetic validation from real-data validation | The manuscript states that biological validity is unestablished and pre-specifies the remaining real-data study |

## Remaining evidence before submission

1. Independently annotated sections from several technologies and tissue types.
2. Multiple upstream clusterers with every post-processor applied to identical
   initial labels.
3. Validation-only tuning and hierarchical inference at subject or slide level.
4. Histology or expert review of changed assignments, especially boundaries.
5. Marker separation, spatially variable genes, cross-section stability, and
   complete reporting of failures.
