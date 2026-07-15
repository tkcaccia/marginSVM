\name{spatial_backend_capabilities}
\alias{spatial_backend_capabilities}
\alias{register_spatial_backend}
\title{Optional Spatial Refinement Backends}
\usage{
spatial_backend_capabilities()

register_spatial_backend(backend = c("cuda", "metal"), provider)
}
\arguments{
\item{backend}{Either \code{"cuda"} or \code{"metal"}.}
\item{provider}{A backend function, or \code{NULL} to unregister it. The function
receives \code{xy}, integer \code{labels}, and a validated \code{control} list.
It returns integer labels with \code{confidence}, \code{changes}, and
\code{neighbors} attributes.}
}
\value{\code{spatial_backend_capabilities()} returns a data frame.
\code{register_spatial_backend()} invisibly returns the backend name.}
\description{
Reports built-in and registered execution backends or registers an optional CUDA
or Metal provider. The package always includes the C++ CPU backend. Accelerator
requests fall back to CPU with a warning unless \code{strict_backend = TRUE}.
}
