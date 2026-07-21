\name{corrupt_categorical_mask}
\alias{corrupt_categorical_mask}
\title{Generate a corrupted categorical mask}
\usage{
corrupt_categorical_mask(
  mask,
  mechanism = c("impulse", "boundary", "patch"),
  rate = 0.15,
  seed = 1L
)
}
\arguments{
\item{mask}{A reference categorical matrix or three-dimensional array.}
\item{mechanism}{One of \code{"impulse"}, \code{"boundary"}, or \code{"patch"}.}
\item{rate}{Fraction of non-missing cells whose labels are replaced.}
\item{seed}{One integer random seed.}
}
\value{A mask with the same representation as \code{mask} and a reproducible set of changed labels. Exactly \code{round(rate * n)} non-missing cells are changed while at least one correct exemplar of every reference class is retained; a rate that makes this impossible is rejected.}
\description{Creates independent, boundary-concentrated, or coherent-patch errors for mask-cleaning experiments.}
