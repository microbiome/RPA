% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/d.update.fast.R
\name{d.update.fast}
\alias{d.update.fast}
\title{Fast d update}
\usage{
d.update.fast(St, s2)
}
\arguments{
\item{St}{probes x samples data matrix}

\item{s2}{variances for the probes}
}
\description{
Computes weighted average over the probes, weighted by their inverse probe-specific variances.
}
\details{
Returns summarized probeset-level weighted average
}
\examples{
# 
}
\references{
See citation("RPA")
}
\author{
Leo Lahti \email{leo.lahti@iki.fi}
}
\keyword{methods}
