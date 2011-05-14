\name{rpa.plot}

\Rdversion{1.1}

\alias{rpa.plot}

\title{Plot RPA results and probe-level data for a specified probeset.}

\description{Plots the preprocessed probe-level observations, estimated
  probeset-level signal, and probe-specific variances. It is also
  possible to highlight individual probes and external summary measures.}

\usage{rpa.plot(dat, rpa.fit.object = NULL, highlight.probes = NULL, pcol =
"darkgrey", mucol = "black", ecol = "red", cex.lab = 1.5, cex.axis = 1, cex.main = 1, cex.names = 1, external.signal = NULL, main = "", plots = "all", ...)
                
}

\arguments{
  \item{dat}{
    Original data: probes x samples.
  }
  \item{rpa.fit.object }{An instance of the 'rpa.fit' class.}

  \item{highlight.probes }{Optionally highlight some of the probes
                            (with dashed line)}
  \item{pcol }{Color for probe signal visualization.}
  \item{mucol }{Color for summary estimate.}
  \item{ecol }{Color for external signal.}  
  \item{cex.lab, cex.axis, cex.main, cex.names}{Font size adjustment parameters.}
  \item{external.signal }{Plot external signal on the probeset. For
                           instance, an alternative summary estimate
                           from another preprocessing methods}
  \item{main}{Title text.}
  \item{plots}{"all": plot data and summary, noise and affinity
               "data": plot data and summary}
  \item{...}{Other parameters to pass for plot function.}
}

\value{Used for its side-effects. Returns probes x samples matrix of
  probe-level data plotted on the image.}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE. See http://www.cis.hut.fi/projects/mi/software/RPA/ }

\author{Leo Lahti \email{leo.lahti@iki.fi}}

\seealso{RPA.pointestimate}

\examples{

# Not run:

## Load example data set
#require(affydata)
#data(Dilution)

## Compute RPA for specific probesets only	
#set <- "1000_at"
#rpa.results <- RPA.pointestimate(Dilution, set)

## Visualize the results for one of the probe sets
#rpa.plot(set, rpa.results)

}

\keyword{ methods }

