\name{plot.rpa}
\alias{plot.rpa}
\title{plot.rpa
Plot RPA results and probe-level data for a specific probeset.}
\usage{
  \method{plot}{rpa}(x, y = NULL, set, highlight.probes = NULL,
    pcol = "darkgrey", mucol = "black", ecol = "red",
    cex.lab = 1.5, cex.axis = 1, cex.names = 1,
    cex.main = 1, external.signal = NULL, main = NULL,
    plots = "all", ...)
}
\arguments{
  \item{x}{Instances of the 'rpa class;}

  \item{y}{Instances of the 'rpa class; y is optional and
  never used. Provided for consistency.}

  \item{set}{Probeset to plot.}

  \item{highlight.probes}{Optionally highlight some of the
  probes (with dashed line)}

  \item{pcol}{Color for probe signal visualization.}

  \item{mucol}{Color for summary estimate.}

  \item{ecol}{Color for external signal.}

  \item{cex.lab}{Label adjustment parameters.}

  \item{cex.axis}{Axis adjustment parameters.}

  \item{cex.names}{Names adjustment parameters.}

  \item{cex.main}{Title adjustment parameters.}

  \item{external.signal}{Plot external signal on the
  probeset. For instance, an alternative summary estimate
  from another preprocessing methods.}

  \item{main}{Title for plot.}

  \item{plots}{Type of plot. See help(rpa.plot).}

  \item{...}{Other arguments to be passed.}
}
\value{
  Used for side-effects. Returns probes x samples matrix of
  probe-level data plotted on the image.
}
\description{
  plot.rpa Plot RPA results and probe-level data for a
  specific probeset.
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
\keyword{utilities}

