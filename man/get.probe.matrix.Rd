\name{get.probe.matrix}
\alias{get.probe.matrix}
\title{get.probe.matrix}
\usage{
  get.probe.matrix(cels, cdf = NULL, quantile.basis,
    bg.method = "rma", normalization.method = "quantiles",
    batch = NULL, verbose = TRUE)
}
\arguments{
  \item{cels}{List of CEL files to preprocess}

  \item{cdf}{Specify an alternative CDF environment}

  \item{quantile.basis}{Pre-calculated basis for quantile
  normalization}

  \item{bg.method}{Specify background correction method.
  See bgcorrect.methods() for options.}

  \item{normalization.method}{normalization method}

  \item{batch}{batch}

  \item{verbose}{Print progress information during
  computation}
}
\description{
  get.probe.matrix
}
\details{
  Returns background-corrected, quantile normalized log2
  probes x samples matrix
}
\examples{
#
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("RPA")
}
\keyword{methods}

