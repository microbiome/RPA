\name{get.probe.parameters}

\Rdversion{1.1}

\alias{get.probe.parameters}

\title{Collect probe-level hyperparameter from batch files}

\description{Collect probe-level hyperparameter from batch files}

\usage{get.probe.parameters(affinities, unique.run.identifier, save.batches.dir = ".", mode = "list")}

\arguments{

  \item{ affinities }{affinities}
  \item{ unique.run.identifier }{Batch file identifier string}
  \item{ save.batches.dir }{Batch file directory}
  \item{ mode  }{Return format: "list" or "table"}

}

\details{Collect parameters from batch files}

\value{List (mode = "list") or data.frame (mode = "table") }

\references{See citation("RPA")}

\author{Leo Lahti \email{leo.lahti@iki.fi}}

\seealso{RPA.pointestimate}

\examples{
# TBA
}

\keyword{ methods }

