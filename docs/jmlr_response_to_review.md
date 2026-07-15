# Response to the JMLR-style review

We thank the reviewer for identifying weaknesses in task formalization, experimental
range, component attribution, and systems evidence. We revised the manuscript and
package rather than treating the new simulations as a cosmetic supplement.

## Major comments

### 1. Learning task and identifiability

**Addressed.** We now define coordinates $X$, observed labels $\tilde y$, latent
reference labels $y^*$, and tissue blocks $g$. The Methods state the piecewise-local
regularity regime and give a non-identifiability argument for coherent component
relabeling. The paper now distinguishes recoverable pointwise/boundary errors from
changes that require expression, histology, or external annotation evidence.

### 2. Novelty, pseudocode, and ablation

**Addressed with a qualified conclusion.** The Introduction now distinguishes
marginSVM from localized SVMs and the earlier KODAMA spatial SVM premise. Methods
include complete pseudocode. We added C++/R ablation switches for cross-fitting and
adaptive tiles and ran 60 held-out scenarios across five structures, two density
profiles, and two error mechanisms. TV decoding and graph fusion made the clearest
contributions; adaptive tiles made a smaller contribution. Cross-fitting and the
robust ramp were nearly neutral, and zero overlap was marginally higher on average.
The Discussion no longer attributes the full gain to every component.

### 3. Corruption fraction

**Addressed.** We ran 792 new datasets spanning all 11 geometries, uniform and
extreme concentration, four error mechanisms, 5%, 15%, and 40% corruption, and
three seeds. Joining the matched 25% block produced 1,056 datasets. marginSVM was
first at 15%, 25%, and 40%, but not at 5%; SpaGCN and graph refinement were safer
when almost all labels were correct. The abstract, Results, and Discussion now make
the operating recommendation conditional on expected error burden.

### 4. Statistical units

**Addressed.** Three seeds are averaged within each of 236 unique simulation
conditions before inference. Paired bootstrap resampling is stratified by geometric
versus gradient family. The paper reports mean rank, win/loss fractions, Friedman
and Holm-adjusted paired tests, and condition-level effect intervals. It explicitly
limits these inferences to the simulator rather than treating the p-values as proof
of biological generalization.

### 5. Rare-region behavior

**Addressed.** We added sparsest-region accuracy as a primary diagnostic, a figure
across all concentration/error combinations, and explicit loss conditions. SpaGCN
achieved 0.7409 sparse-region accuracy versus 0.6644 for marginSVM. The Discussion
now recommends conservative refinement when very rare compartments dominate the
scientific objective.

### 6. Hyperparameters and trust-field mode

**Addressed in part.** The default and experimental modes are separated throughout.
The trust-field extension remains opt-in because it improved the frozen 19-class
colorectal block but lost to the default in the 708-scenario matrix. Component
ablation and the existing disjoint development/confirmation split are reported.
A full high-dimensional sensitivity surface would be expensive and is not claimed.

### 7. Comparator validity

**Addressed.** The manuscript formalizes input equivalence: direct comparators must
consume only coordinates, labels, and tissue identifiers. End-to-end clustering is
excluded because it consumes expression or image data and estimates labels from
scratch. Package tests compare the native SpaGCN and GraphST rules against literal
small-data reference implementations, including neighborhood sizes and tie behavior.

### 8. Complexity and systems evidence

**Addressed.** Methods now give time and memory complexity in $n$, tile halo sizes,
landmarks, classes, graph degree, and TV iterations. A fresh marginSVM scaling study
uses separate processes and measures peak resident memory. At 500,000 observations,
runtime was 10.28 seconds/752 MB in 2D and 15.48 seconds/841 MB in 3D with four CPU
workers. CUDA and Metal are described only as provider contracts; no accelerator
performance claim remains.

### 9. Independent biological validation

**Outstanding by design.** The colorectal experiment recovers a reference
annotation after synthetic corruption; it does not establish biological truth. We
have narrowed the paper to algorithm and software validation and state that DLPFC,
additional tumors, platforms, marker separation, histological concordance, and
expert review remain necessary before a biological methods claim. This limitation
cannot be removed without the independent datasets.

## Minor comments

The manuscript now uses `marginSVM` consistently, calls softmax outputs evidence
rather than calibrated probabilities, labels 3D panels as projected displays,
clarifies tissue isolation, lists exact regeneration commands, and removes
request-oriented or externally supplied-data wording. Raw scenario-level CSV files,
seeds, package tests, and all figure scripts are retained with the source.
