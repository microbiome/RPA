\name{rpa}

\Rdversion{1.1}

\alias{rpa}
\alias{rpa.online}

\title{RPA for preprocessing.}

\description{Returns an expressionSet object preprocessed with RPA. If
  'cind' is not specified, uses the first array of affybatch as the
  reference. With rpa.online one can run preprocessing in online-mode
  to avoid memory problems with huge expression data collections.}

\usage{

rpa(abatch = NULL, sets = NULL, priors = NULL, epsilon = 1e-2, cind = 1, sigma2.method = "robust", d.method = "fast", verbose = FALSE, bg.method = "rma", normalization.method = "quantiles.robust", cdf = NULL, cel.files = NULL, cel.path = NULL) 

rpa.online(cel.path = NULL, cel.files = NULL, sets = NULL, priors = list(alpha = 2, beta = 1), epsilon = 0.01, cind = 1, verbose = FALSE, bg.method = "rma", cdf = NULL, batch.size = 10)

}

\arguments{

  \item{abatch }{An AffyBatch object.}

  \item{sets }{
    Probesets for which RPA will be computed. Default: all probe sets.}

  \item{priors }{An 'rpa.priors' object. Can be used to set
		user-specified priors for the model parameters. Not
		used sigma2.method = "var". The prior parameters alpha
		and beta are prior parameters for inverse Gamma
		distribution of probe-specific
		variances. Noninformative prior is obtained with
		alpha, beta -> 0.  Not used with sigma2.method
		'var'. Scalar alpha and beta specify an identical
		inverse Gamma prior for all probes, which regularizes
		the solution. Can be also specified as lists, each
		element corresponding to one probeset.}

  \item{epsilon }{Convergence tolerance. The iteration is deemed converged 
		  when the change in all parameters is < epsilon.}
  
  \item{cind }{Specify reference array for computing probe-level
    	       differential expression. Default: cind = 1. Note that if
    	       exclude.reference.array = TRUE the expression value for
    	       the reference array (cind) will be excluded in the
    	       output. Note that all values of the reference array are 0
    	       since they indicate the differential expression of the
    	       reference array against itself.}

  \item{sigma2.method }{

    Optimization method for sigma2 (probe-specific variances). This
    parameter is denoted by tau^2 in the vignette and manuscript.

	"robust": (default) update sigma2 by posterior mean,
		regularized by informative priors that are identical
		for all probes (user-specified by
		setting scalar values for alpha, beta). This
		regularizes the solution, and avoids overfitting where
		a single probe obtains infinite reliability. This is a
	        potential problem in the other sigma2 update
	        methods with non-informative variance priors. The
		default values alpha = 2; beta = 1 are
	        used if alpha and beta are not specified.
   
        "mode": update sigma2 with posterior mean

	"mean": update sigma2 with posterior mean
	
	"var": update sigma2 with variance around d. Applies the fact
               that sigma2 cost function converges to variance with
               large sample sizes. 

		}

  \item{d.method }{
    Method to optimize d.


        "fast": (default) weighted mean over the probes, weighted by
		probe variances The solution converges to this with
		large sample size.

        "basic": optimization scheme to find a mode used in Lahti et
        	 al. TCBB/IEEE; relatively slow; this is the preferred 
		 method with small sample sizes.
                
	      }

    \item{verbose }{Print progress information during computation.}

    \item{bg.method }{
      Specify background correction method. Default:
      "rma". See bgcorrect.methods() for other options.}
  
    \item{normalization.method }{
      Specify quantile normalization method. Default: "pmonly". See
      normalize.methods(Dilution) for other options.}

    \item{cdf }{
      Specify an alternative CDF environment. Default: none.
    }

    \item{cel.files}{List of CEL files to preprocess.}

    \item{cel.path}{Path to CEL file directory.}

    \item{batch.size}{Batch size for online mode (rpa.online); the complete list of CEL files will be preprocessed in batches with this size using Bayesian online-updates for probe-specific parameters.}
}

\details{RPA preprocessing function. Gives an estimate of the
  probeset-level mean parameter d of the RPA model, and returns these
  in an expressionSet object. In online-mode (rpa.online), the CEL
  files are handled in batches to obtain Bayesian updates for
  probe-specific hyperpriors; after scanning through the database the
  results are combined. The online mode is useful for preprocessing
  very large expression data sets where ordinary preprocessing
  algorithms fail.}

\value{
  An instance of the 'expressionSet' class.
}

\references{See citation("RPA").}

\author{Leo Lahti \email{leo.lahti@iki.fi}}

\note{sigma2.method = "robust" and d.method = "fast" are
recommended. With small sample size and informative prior, d.method =
"basic" may be preferable.}

\seealso{RPA.pointestimate, set.priors, AffyBatch, ExpressionSet,
  estimate.affinities, rpa.fit}

\examples{

# Not run:

## Load example data set
#require(affydata)
#data(Dilution)

## Compute RPA for specific probesets
#sets <- geneNames(Dilution)[1:2]	
#set <- "33572_at"
#eset <- rpa(Dilution, sets)

## Compute RPA for whole data set	
## ... slow, not executed here	
## eset <- rpa(Dilution)

}

\keyword{ methods }

