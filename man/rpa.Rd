\name{rpa}

\Rdversion{1.1}

\alias{rpa}

\title{RPA for preprocessing.}

\description{Returns an expressionSet object preprocessed with RPA. If
  'cind' is not specified, uses the first array of affybatch as the
  reference.}

\usage{rpa( abatch,
                    sets = NULL,
                  myseed = 101, 
                  priors = NULL,
                 epsilon = 1e-2, 
                    cind = 1,
           sigma2.method = "robust",
                d.method = "fast",
                 verbose = FALSE,
               bg.method = "rma",
    normalization.method = "quantiles.robust",
                     cdf = NULL,
                   alpha = NULL,
                    beta = NULL,
         affinity.method = "rpa", quantile.n = 50) 
}

\arguments{

  \item{abatch }{An AffyBatch object.}

  \item{sets }{
    Probesets for which RPA will be computed. Default: all probe sets.}

  \item{myseed }{Specify random seed.}

  \item{priors }{An 'rpa.priors' object. Can be used to set
		user-specified priors for the model parameters. Not used
		sigma2.method = "var".}

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

    \item{alpha, beta }{Prior (scalar) parameters for inverse Gamma
      distribution of probe-specific variances. Noninformative prior is
      obtained with alpha, beta -> 0.  Not used with sigma2.method
      'var'. Scalar alpha and beta specify an identical inverse Gamma
      prior for all probes, which regularizes the
      solution. Probe-specific priors can be set with the 'priors'
      parameter.}

    %\item{exclude.reference.array }{Logical indicating whether the
    %  values for the reference array will be excluded in the final
    %  differential gene expression matrix. Note that all values of the
    %  reference array will be 0 since they indicate differential
    %  expression of the reference array against itself.}

    \item{affinity.method}{
      For model details, see 'help(estimate.affinities)'.
      
      "rpa": Assuming affinity parameters are zero on average, and the
      deviation from zero is determined by estimated probe-level noise
      parameters. This gives higher weight (smaller affinity) for more
      reliable probes also in affinity estimation. Heuristic solution,
      which aims to fit probe-level signal in real data domain as close
      to the reliable probes as possible.

      "zeromean": assumes that probe affinities sum to zero. Analogous
      to model assumptions in RMA. Gives equal weights for all probes in
      affinity estimation. We expect this to be less optimal than
      weighting probes by their general reliability.
      
    }

    \item{quantile.n }{Number of samples to be used with quantiles.online normalization method.}

}

\details{RPA preprocessing function. Gives an estimate of the
  probeset-level mean parameter d of the RPA model, and returns these in
  an expressionSet object.}

\value{
  An instance of the 'expressionSet' class.
}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE. See
http://www.cis.hut.fi/projects/mi/software/RPA/ }

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

