\name{simulate_spatial_clusters}
\alias{simulate_spatial_clusters}
\title{Simulate labelled spatial clusters}
\usage{
simulate_spatial_clusters(
  n = 50000L,
  dimensions = 2L,
  k = 5L,
  samples = 1L,
  noise = 0.08,
  seed = 1L
)
}
\arguments{
\item{n}{Number of observations.}
\item{dimensions}{Spatial dimensions, usually 2 or 3.}
\item{k}{Number of clusters.}
\item{samples}{Number of independent samples or slides.}
\item{noise}{Fraction of labels to randomly corrupt.}
\item{seed}{Random seed.}
}
\value{
A list with \code{xy}, \code{labels}, \code{truth}, and \code{samples}.
}
\description{
Creates a simple spatial benchmark dataset with separated domains and optional
label corruption.
}
\examples{
sim <- simulate_spatial_clusters(n = 500, dimensions = 3, k = 4)
str(sim)
}
