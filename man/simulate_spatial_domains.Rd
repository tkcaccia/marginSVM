\name{simulate_spatial_domains}
\alias{simulate_spatial_domains}
\title{Simulate Spatial Tissue Domains with Challenging Geometry}
\usage{
simulate_spatial_domains(
  n = 50000L,
  pattern = c("jagged_stripes", "wavy_layers", "rings", "spiral",
    "branching", "lobes", "islands", "disconnected", "thin_layers",
    "intermixed", "layers3d"),
  k = 5L,
  noise = 0.2,
  dimensions = if (identical(pattern[1L], "layers3d")) 3L else 2L,
  samples = 1L,
  noise_type = c("random", "boundary", "patch", "region"),
  feature_scale = 1,
  density_profile = c("uniform", "moderate", "strong", "extreme", "hotspot"),
  seed = 1L
)
}
\arguments{
\item{n}{Number of observations.}
\item{pattern}{Tissue-domain geometry.}
\item{k}{Number of domains.}
\item{noise}{Fraction of initial labels to corrupt.}
\item{dimensions}{Two or three spatial dimensions.}
\item{samples}{Number of independent samples.}
\item{noise_type}{One of \code{"random"}, \code{"boundary"},
\code{"patch"}, or \code{"region"}.}
\item{feature_scale}{Relative width of islands and thin layers.}
\item{density_profile}{Relative observation concentration across tissue
regions. Use \code{"uniform"}, \code{"moderate"}, \code{"strong"},
\code{"extreme"}, \code{"hotspot"}, or a positive numeric vector with one
weight per region. Character profiles reproducibly permute the weights across
region identifiers.}
\item{seed}{Random seed.}
}
\value{A list with coordinates, noisy labels, truth, region and sample
identifiers, realized region counts, concentration weights, and simulation
metadata.}
\description{
Generates irregular tissue sections or volumes containing curved, thin, rare,
or disconnected spatial domains for evaluating cluster-label refinement.
}
