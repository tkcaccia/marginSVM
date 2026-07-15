# Literature And Tool Context For Spatial Domain Benchmarks

The benchmark suite in `benchmarks/benchmark_suite.R` follows common dimensions
used by recent spatial transcriptomics method evaluations:

- Accuracy against ground truth annotations, usually reported with ARI or related
  clustering agreement metrics.
- Spatial continuity or smoothness of predicted domains.
- Robustness to noise, platform/resolution, and biological replicate/sample
  structure.
- Scalability with increasing spot/cell counts.
- Visual inspection of spatial maps, because similar numeric scores can still
  produce visibly different tissue boundaries.

End-to-end spatial domain tools reported in recent benchmarks include BayesSpace,
BASS, BANKSY, PRECAST, SpaGCN, STAGATE, GraphST, SEDR, DeepST, DR.SC,
SpatialPCA, Leiden, Louvain, and Walktrap. They are useful context but are not
comparators for this coordinate-and-label post-processing study.

Useful benchmark references:

- Li et al., Nature Methods 2024: benchmarked 13 methods on 34 SRT sections and
  evaluated accuracy, spatial continuity, marker genes, scalability, and
  robustness. DOI: 10.1038/s41592-024-02215-8
- Kang et al., Nucleic Acids Research 2025: benchmarked 19 methods across 30
  real datasets and 27 synthetic datasets, focusing on domain detection and
  domain-specific spatially variable genes. DOI: 10.1093/nar/gkaf303
- SRTBenchmark: benchmarked 14 spatial clustering methods across roughly 600
  datasets, ten technologies, and eight organs.
  https://github.com/ZJUFanLab/SRTBenchmark
This package is currently benchmarked as a post-processing/refinement method,
not as a full end-to-end domain discovery tool. Direct competitors are the
published SpaGCN strict-majority and GraphST nearest-neighbor modal refinement
functions, matched neighborhood voting, Potts-like smoothing, and morphology on
regular grids. The precursor KODAMA SVM and all end-to-end clustering pipelines
are deliberately excluded from comparative claims.
