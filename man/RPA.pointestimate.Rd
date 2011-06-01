\name{RPA.pointestimate}

\Rdversion{1.1}

\alias{RPA.pointestimate}

\title{Computing point estimate for the model parameters for all probe
sets.}

\description{Computes point estimate}

\usage{RPA.pointestimate(abatch, sets = NULL, myseed = 101, priors =
NULL, epsilon = 1e-2, cind = 1, sigma2.method = "robust", d.method =
"fast", verbose = TRUE, bg.method = "rma", normalization.method =
"quantiles.robust", cdf = NULL, alpha = NULL, beta = NULL, quantile.n = 50)}

\arguments{

  \item{abatch }{An AffyBatch object.}

  \item{sets }{Specifies the probesets for which RPA estimates will be
  computed. Default: all probe sets.}

  \item{myseed }{Specifies the random seed.}

  \item{priors }{An 'rpa.priors' object. Can be used to set
		 user-specified priors for the model parameters. Not
		 applicable for sigma2.method = "var".}
  
  \item{epsilon }{Convergence tolerance. The iteration is deemed
  	 	  converged when the change in the d parameter 
		  is < epsilon.}
  
  \item{cind }{Specifies which array in abatch is used as a
	       reference in computing probe-level differential
	       expression.}

  \item{sigma2.method }{

    Optimization method for sigma2 (probe-specific variances).

	"robust": (default) update sigma2 by posterior mean,
		regularized by informative priors that are identical
		for all probes (user-specified by
		setting scalar values for alpha, beta). This
		regularizes the solution and avoids overfitting where
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

    \item{verbose }{
      Print progress information during computation. Default: TRUE.}

    \item{bg.method }{
      Specify background correction method. Default:
      "rma". See bgcorrect.methods() for other options.}
  
    \item{normalization.method }{
      Specify quantile normalization method. Default: "pmonly". See
      normalize.methods(Dilution) for other options.}

    \item{cdf }{
      Specify an alternative CDF environment. Default: none.
    }

    \item{alpha, beta }{Prior parameters for inverse Gamma
      distribution of probe-specific variances. Noninformative prior is
      obtained with alpha, beta -> 0.  Not used with sigma2.method
      'var'. Scalar alpha and beta specify an identical inverse Gamma
      prior for all probes, which regularizes the
      solution. Probe-specific priors can be set with the 'priors'
      parameter.}

    \item{quantile.n }{Number of samples to be used with quantiles.online normalization method.}

}

\details{Calculates RPA estimates of probe reliability and differential
  expression between the user-specified reference array (cind) and the
  other arrays in the data set. The model assumes P observations for
  each transcript target (i.e. a probeset) with Gaussian noise which is
  specific for each probe (variance is specified by sigma2). The mean
  (affinity) parameters of the Gaussian noise model cancel out in
  calculating probe-level differential expression.  RPA.pointestimate
  gives a point estimate for d and sigma2. The 'prior' parameter is not
  applicable with sigma2.method = "var". The d.method = "fast" is
  recommended with large sample size.}

\value{

  An instance of class 'rpa'. This is an extended list containing
  the following elements:

  \item{d }{A matrix of probesets x arrays. Specifies the estimated
            'true' underlying differential gene expression signal over
	    the arrays (vs. the reference array 'cind') for each 
	    investigated probeset. Note that the reference array is 
	    not included.}
  
  \item{sigma2 }{A list. Each element corresponds to a probeset, and
                 contains a vector that gives the estimated variance for
                 each probe in that probeset. This corresponds to the
                 parameter tau^2 in the vignette and manuscript.}
  
  \item{cind }{Specifies which of the arrays in the abatch (the
	       affybatch object to be analyzed) has been used as the 
	       reference for computing probe-level differential
	       expression.}

  \item{affinity}{Probe affinity effects.}
	     
  \item{sets }{A character vector listing the investigated probesets.}
}

\references{See citation("RPA").}

\author{Leo Lahti \email{leo.lahti@iki.fi}}

\note{sigma2.method = "robust" and d.method = "fast" are
recommended. With small sample size and informative priors, d.method =
"basic" may be preferable, with large sample size d.method = "fast"
should have considerable speedups and comparable accuracy.}

\seealso{rpa.plot, rpa, set.priors, rpa2eset, RPA.preprocess, AffyBatch,
rpa.fit, estimate.affinities}

\examples{

## Load example data set
#require(affydata)
#data(Dilution)

## Compute RPA for whole data set	
## ... slow, not executed here	
# rpa.results <- RPA.pointestimate(Dilution)

## Visualize the results for one of the probe sets
#rpa.plot(set, rpa.results)

}

\keyword{ methods }

