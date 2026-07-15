# marginSVM v2: Local Trust and Boundary Fields

Working name: marginSVM.

## Implementation Status (July 2026)

The five proposed mechanisms are implemented in the C++17 engine behind the
opt-in `control = list(experimental_v2 = 1)` path. The validated v1 path remains
the default. The implementation uses 6/12/24/48-neighbor support, cross-fitted
input-label probability, overlapping-tile variance, and local probability-field
stability to estimate a continuous trust field. A local regularized covariance
reranks kd-tree candidates with a Euclidean/Mahalanobis mixture. Boundary points
receive bounded one-vs-one Nyström SVM adjustments for the most frequent ambiguous
class pairs. Output attributes report trust, tile disagreement, perturbation
stability, selective risk, retain/change/unresolved state, and protected rare
components.

Development on eight 194,541-point, 19-class colorectal corruptions used seed
block 840000. With all controls frozen, a disjoint confirmation block beginning at
1040000 improved mean accuracy from 0.8677 (v1) to 0.8717, ARI from 0.7923 to
0.7997, and worst-class recall from 0.5874 to 0.6162 while reducing damage from
0.0232 to 0.0128. Mean boundary accuracy changed from 0.7383 to 0.7361 and runtime
increased from 4.61 to 7.26 seconds.

The extension did not advance as a universal replacement: on 16 mixed 2D/3D
development simulations, v2 accuracy was 0.9005 versus 0.9224 for v1. It therefore
remains an explicitly experimental high-class mode. This negative result prevents
selection on the colorectal tissue from being generalized to all geometries.

## Motivation

The current method is strongest overall but still joins two behaviors with a
high-cardinality rule: aggressive correction of dispersed errors and tissue-level
abstention for coherent changes. A better model should infer correction strength
continuously at each observation and should distinguish a real thin boundary from
random label discordance.

## Proposed Model

For observation `i`, estimate two latent fields:

- `r_i`: probability that the input label is reliable;
- `b_i`: probability that the observation lies on a biological boundary.

The reliability feature vector contains only out-of-fold quantities:

1. cross-fitted SVM probability of the input label;
2. same-label support at 6, 12, 24, and 48 neighbors;
3. entropy and Jensen-Shannon disagreement of neighboring SVM probabilities;
4. variance across regular/adaptive tiles and two landmark resamples;
5. class prevalence and local connected-component size;
6. prediction stability under small coordinate perturbations.

Fit a constrained two-component mixture to these features without reference
truth. The reliable component must have higher multiscale support and higher
perturbation stability; class-stratified priors prevent rare classes from being
assigned wholesale to the unreliable component.

The boundary field combines the spatial gradient of cross-fitted scores, local
probability divergence, and anisotropic neighborhood geometry. A local covariance
matrix defines a Mahalanobis distance so edges follow tissue fibers and thin layers
instead of crossing them isotropically.

The joint objective becomes

```
sum_i r_i * capped_hinge_i
+ sum_i,c [-log(P_ic) + logit(r_i) * I(c != y_i)] * u_ic
+ beta * sum_(i,j) (1 - b_ij) * w_ij * ||u_i - u_j||_1
+ eta * sum_i KL(r_i || r_i_previous).
```

High reliability protects the input assignment. High boundary probability weakens
spatial coupling but does not automatically trigger a label change. The KL term
prevents an iterative fit from making abrupt self-confirming reliability changes.

## Boundary Specialists

Train multiclass SVMs only for broad candidate generation. Within the uncertain
boundary belt, fit small one-vs-one Nyström experts for the two or three locally
plausible classes. This should improve the current weakness under boundary-targeted
errors while reducing work in stable interiors.

Pairwise votes should be coupled into probabilities before TV decoding. A boundary
expert is discarded when either class has insufficient reliable support in its
training halo.

## Selective Prediction

Calibrate out-of-fold margins by tissue and class using temperature scaling. Return
three states internally: retain, change, or unresolved. The public factor retains
the input label for unresolved points, while attributes report selective coverage
and estimated risk. This makes abstention pointwise rather than tissue-wide.

## Failure Safeguards

- Never train and score the same observation in one fold.
- Limit reliability/SVM alternation to two iterations.
- Require the structured objective to decrease after an iteration.
- Preserve a minimum reliable mass per initial class.
- Do not merge connected components solely because they are small.
- Freeze all v2 settings before opening a new confirmatory seed block.

## Development Experiments

Use new development seeds only. Compare:

1. current marginSVM;
2. continuous trust field;
3. trust plus anisotropic graph;
4. trust plus pairwise boundary specialists;
5. complete v2 model.

Primary objective: mean accuracy plus 0.3 macro recall plus 0.15 worst-class recall,
with a hard damage-rate ceiling. Report boundary accuracy, correction recall,
selective-risk curves, runtime, peak memory, and per-class connected-component
retention.

The colorectal criterion was met on a fresh block, but the mixed-geometry criterion
was not. Consequently v2 remains opt-in and requires broader independent data
before promotion to the package default.
