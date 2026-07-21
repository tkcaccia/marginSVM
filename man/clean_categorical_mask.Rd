\name{clean_categorical_mask}
\alias{clean_categorical_mask}
\title{Clean a noisy categorical mask}
\usage{
clean_categorical_mask(mask, samples = NULL, workers = NULL)
}
\arguments{
\item{mask}{A categorical matrix or three-dimensional array. Missing cells are treated as void locations.}
\item{samples}{Optional matrix or array with the same dimensions as \code{mask}. Different samples are cleaned independently.}
\item{workers}{Optional CPU budget used across samples and independent spatial charts. \code{NULL} uses up to four physical cores when possible.}
}
\value{A cleaned mask with the representation of \code{mask}. Pointwise
attributes are \code{candidate}, \code{margin_score}, \code{required},
\code{repair_margin}, \code{atlas_dispersion}, \code{isolation}, and
\code{changed}. \code{margin_score} is a local support contrast on the mask
lattice; \code{repair_margin} is the gap to the local acceptance threshold.}
\description{Applies FiberMargin to a 2D pixel mask or 3D voxel mask using only coordinates and categorical assignments.}
\examples{
mask <- matrix(c(rep("A", 50), rep("B", 50)), nrow = 10)
mask[5, 5] <- "B"
cleaned <- clean_categorical_mask(mask, workers = 1)
}
