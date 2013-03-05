\name{rpa.summarize}
\alias{rpa.summarize}
\title{rpa.summarize}
\usage{
  rpa.summarize(dat, affinities, variances,
    summarize.with.affinities = FALSE)
}
\arguments{
  \item{dat}{Original data: probes x samples.}

  \item{affinities}{Probe affinities}

  \item{variances}{Probe variances}

  \item{summarize.with.affinities}{Use affinity estimates
  in probe summarization step. Default: FALSE.}
}
\value{
  A vector. Probeset-level summary signal.
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
# res <- rpa.summarize(dat, affinities, variances, summarize.with.affinities = FALSE)
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("RPA")
}
\seealso{
  rpa
}
\keyword{utilities}

