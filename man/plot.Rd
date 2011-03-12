\name{plot}

\Rdversion{1.1}

\alias{plot.rpa}

\title{Plot RPA results and probe-level data for a specific probeset.}

\description{Plots the preprocessed probe-level observations, estimated
  probeset-level signal, and probe-specific variances. It is also
  possible to highlight individual probes and external summary measures.}

\usage{plot.rpa(x, y = NULL, set, highlight.probes = NULL, pcol =
"darkgrey", mucol = "black", ecol = "red", cex.lab = 1.5, cex.axis = 1,
external.signal = NULL, main = NULL, ...)}

\arguments{
  \item{x,y }{Instances of the 'rpa class; y is optional and never
  used. Provided for consistency.}
  \item{set}{Probeset to plot.}
  \item{highlight.probes }{Optionally highlight some of the probes
                           (with dashed line)}
  \item{pcol }{Color for probe signal visualization.}
  \item{mucol }{Color for summary estimate.}
  \item{ecol }{Color for external signal.}  
  \item{cex.lab, cex.axis }{Axis adjustment parameters.}
  \item{external.signal }{Plot external signal on the probeset. For
                           instance, an alternative summary estimate
                           from another preprocessing methods.}
  \item{main}{Title for plot.}			 
  \item{...}{Other arguments to be passed.}
			 
}		       

\value{Used for its side-effects. Returns probes x samples matrix of
  probe-level data plotted on the image.}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE. See http://bioconductor.org/packages/release/bioc/html/RPA.html}

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
#plot.rpa(set, rpa.results)

}

\keyword{ methods }

