\name{refine_spatial_clusters}
\alias{refine_spatial_clusters}
\title{Refine Spatial Cluster Assignments}
\usage{
refine_spatial_clusters(xy, labels, samples = NULL, control = NULL, execution = NULL)
}
\arguments{
\item{xy}{Numeric matrix of two- or three-dimensional spatial coordinates.}
\item{labels}{Initial cluster assignment for each observation.}
\item{samples}{Optional sample or slide identifier.}
\item{execution}{Optional named list controlling large-data execution. Use
\code{tiles = "auto"} for automatically sized overlapping tiles. Other entries
are \code{overlap}, \code{target_tile_size}, \code{workers}, \code{backend}, and
\code{strict_backend}.}
\item{control}{Optional named list with advanced settings: \code{neighbors},
\code{iterations}, \code{consensus}, \code{preserve}, \code{margin}, and
\code{weighted}, and \code{current_support}.}
}
\value{A factor of refined assignments with local support and change diagnostics.
The support value is not a calibrated probability.}
\description{
Builds an exact nearest-neighbor graph in C++ and synchronously repairs locally
inconsistent labels while retaining assignments at ambiguous tissue boundaries.
Defaults are inferred from dataset size. Tile cores do not overlap, while their
input halos do; only core predictions are retained so each observation has one
deterministic owner.
}
\examples{
sim <- simulate_spatial_domains(2000, "jagged_stripes")
refined <- refine_spatial_clusters(sim$xy, sim$labels)
mean(refined == sim$truth)

tiled <- refine_spatial_clusters(
  sim$xy, sim$labels,
  execution = list(tiles = c(2, 2), workers = 2)
)
}
