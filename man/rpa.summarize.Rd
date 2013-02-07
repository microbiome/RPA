\name{rpa.summarize}
\alias{rpa.summarize}
\title{rpa.summarize}
\usage{
  rpa.summarize(dat, affinities, variances)
}
\arguments{
  \item{dat}{Original data: probes x samples.}

  \item{affinities}{Probe affinities}

  \item{variances}{Probe variances}
}
\description{
  rpa.summarize
}
\details{
  Summarizes the probes in a probe set according to the RPA
  model based on the given affinity and variance
  parameters.
}
\examples{
# res <- rpa.summarize(dat, affinities, variances)
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("RPA")
}
\seealso{
  rpa, RPA.pointestimate
}
\keyword{utilities}

