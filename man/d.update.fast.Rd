\name{d.update.fast}
\alias{d.update.fast}
\title{d.update.fast}
\usage{
  d.update.fast(St, s2)
}
\arguments{
  \item{St}{probes x samples data matrix}

  \item{s2}{variances for the probes}
}
\description{
  Description: computes weighted average over the probes,
  weighted by their inverse probe-specific variances.
}
\details{
  Returns summarized probeset-level weighted average
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

