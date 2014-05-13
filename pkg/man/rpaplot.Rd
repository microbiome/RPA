\name{rpaplot}
\alias{rpaplot}
\title{rpaplot
Plot RPA results and probe-level data for a specified probeset.}
\usage{
  rpaplot(dat, mu = NULL, tau2 = NULL, affinity = NULL,
    highlight.probes = NULL, pcol = "darkgrey",
    mucol = "black", ecol = "red", cex.lab = 1.5,
    cex.axis = 1, cex.main = 1, cex.names = 1,
    external.signal = NULL, main = "", plots = "all", ...)
}
\arguments{
  \item{dat}{Background-corrected and normalized data:
  probes x samples.}

  \item{mu}{probeset signal}

  \item{tau2}{probe variances}

  \item{affinity}{probe affinities}

  \item{highlight.probes}{Optionally highlight some of the
  probes (with dashed line)}

  \item{pcol}{Color for probe signal visualization.}

  \item{mucol}{Color for summary estimate.}

  \item{ecol}{Color for external signal.}

  \item{cex.lab}{Label size adjustment parameters.}

  \item{cex.axis}{Axis size adjustment parameters.}

  \item{cex.main}{Title size adjustment parameters.}

  \item{cex.names}{Names size adjustment parameters.}

  \item{external.signal}{Plot external signal on the
  probeset. For instance, an alternative summary estimate
  from another preprocessing methods}

  \item{main}{Title text.}

  \item{plots}{"all": plot data and summary, noise and
  affinity; "data": plot data and summary}

  \item{...}{Other parameters to pass for plot function.}
}
\value{
  Used for its side-effects. Returns probes x samples
  matrix of probe-level data plotted on the image.
}
\description{
  rpaplot Plot RPA results and probe-level data for a
  specified probeset.
}
\details{
  Plots the preprocessed probe-level observations,
  estimated probeset-level signal, and probe-specific
  variances. It is also possible to highlight individual
  probes and external summary measures.
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

