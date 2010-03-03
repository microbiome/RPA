\name{rpa.plot}

\Rdversion{1.1}

\alias{rpa.plot}

\title{Plot RPA results and probe-level data for a specified probeset.}

\description{Plots the preprocessed probe-level observations, estimated
  probeset-level signal d, and probe-specific variances. It is also
  possible to highlight individual probes.}

\usage{rpa.plot(set, rpa.object, highlight.probes = NULL, pcol = "darkgrey", dcol = "black", cex.lab = 1.5, cex.axis = 1)
                
}

\arguments{

  \item{set }{Probeset to visualize.}
    
  \item{rpa.object }{An instance of the 'rpa' class. Provided by
                     'RPA.pointestimate' function.}

  \item{highlight.probes }{Optionally highlight some of the probes
                           (with dashed line)}
  \item{pcol }{Color for probe signal visualization.}
  \item{dcol }{Color for probeset-level summary d.}
  \item{cex.lab, cex.axis }{Axis adjustment parameters.}
}

\value{Used for its side-effects.}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE. See http://www.cis.hut.fi/projects/mi/software/RPA/ }

\author{Leo Lahti <leo.lahti@iki.fi>}

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

