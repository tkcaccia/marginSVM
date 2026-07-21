\name{evaluate_mask_cleaning}
\alias{evaluate_mask_cleaning}
\title{Evaluate categorical mask cleaning}
\usage{
evaluate_mask_cleaning(
  reference,
  initial,
  cleaned,
  elapsed = NA_real_,
  method = NULL
)
}
\arguments{
\item{reference}{Reference categorical matrix or three-dimensional array.}
\item{initial}{Imperfect input mask.}
\item{cleaned}{Output mask to evaluate.}
\item{elapsed}{Optional elapsed runtime in seconds.}
\item{method}{Optional method name included in the returned row.}
}
\value{A one-row data frame containing class-aware overlap, boundary agreement, rare-class, correction, and damage measures.}
\description{Evaluates a cleaned mask against a fixed reference and noisy input mask. The output includes overlap scores, class-boundary agreement, and correction/damage decomposition.}
