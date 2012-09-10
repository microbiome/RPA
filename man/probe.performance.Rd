\name{probe.performance}
\alias{probe.performance}
\title{probe.performance}
\usage{
  probe.performance(rpa.object, abatch = NULL, sets = NULL)
}
\arguments{
  \item{rpa.object}{Object of the rpa class (output from
  function RPA.pointestimate)}

  \item{abatch}{Optional: Affybatch for the rpa object, if
  not provided in the rpa project.}

  \item{sets}{Optional: Specify the inspected probesets.}

}
\value{
  Data frame of probe-level parameter estimates
}
\description{
  Provide a table of probe-level parameter estimates
  (affinity and stochastic noise) for rpa class object.
}
\examples{
# library(affydata) data(Dilution); rpa.results <- RPA.pointestimate(Dilution); tab <- probe.parameters(rpa.results)
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("RPA")
}
\keyword{utilities}

