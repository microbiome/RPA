\name{summarize.batch}
\alias{summarize.batch}
\title{summarize.batch}
\usage{
  summarize.batch(q, set.inds, priors = NULL,
    probe.parameters = NULL, cind = 1, epsilon,
    verbose = FALSE, mc.cores = 1)
}
\arguments{
  \item{q}{Background corrected, quantile-normalized, log2
  probes x samples matrix}

  \item{set.inds}{Indices for each probeset, corresponding
  to q matrix}

  \item{priors}{Optional probe-specific priors}

  \item{probe.parameters}{A list, each element
  corresponding to a probe set. Each probeset element has
  the following elements: affinity, variance. Each of these
  two elements contains a vector over the probeset probes,
  specifying the affinities and variances for the probes
  according to the RPA model. If given, overrides the
  priors.}

  \item{cind}{Specify reference array for computing
  probe-level differential expression. Default: cind = 1.
  Note that if exclude.reference.array = TRUE the
  expression value for the reference array (cind) will be
  excluded in the output. Note that all values of the
  reference array are 0 since they indicate the
  differential expression of the reference array against
  itself.}

  \item{epsilon}{Convergence tolerance. The iteration is
  deemed converged when the change in all parameters is <
  epsilon.}

  \item{verbose}{Print progress information during
  computation.}

  \item{mc.cores}{Number of cores for parallel processing}
}
\description{
  summarize.batch
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

