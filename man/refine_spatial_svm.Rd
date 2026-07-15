\name{refine_spatial_svm}
\alias{refine_spatial_svm}
\title{Refine spatial cluster assignments with marginSVM}
\usage{
refine_spatial_svm(
  xy,
  labels,
  samples = NULL,
  control = NULL,
  backend = c("auto", "cpu", "cuda", "metal"),
  verbose = FALSE
)
}
\arguments{
\item{xy}{Numeric matrix with two or three spatial coordinates per row.}
\item{labels}{Initial cluster assignments.}
\item{samples}{Optional tissue or section identifier per row.}
\item{control}{Optional named list of advanced controls. Most analyses should use the automatic defaults.}
\item{backend}{Computational backend. CPU is always available; accelerator requests fall back explicitly when unavailable.}
\item{verbose}{Print native tile progress.}
}
\value{
A factor with confidence, margin, local-support, tile, backend, worker, and
abstention diagnostics stored as attributes. The experimental v2 path also
returns trust, tile disagreement, perturbation stability, selective risk,
pointwise decision, and protected-component attributes.
}
\description{
Fits robust nonlinear SVM boundaries in overlapping regular and adaptive tiles,
then reconciles cross-fitted scores with uncertainty-gated graph evidence
and edge-aware spatial total variation. Neighborhoods, tiling, Nyström features,
SVM fitting, fusion, and decoding execute in C++; this R function validates inputs
and selects automatic defaults.
}
\examples{
sim <- simulate_spatial_clusters(n = 1000, dimensions = 2)
refined <- refine_spatial_svm(sim$xy, sim$labels)
table(refined)
}
