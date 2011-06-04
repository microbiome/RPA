\name{RPA.iteration}

\Rdversion{1.1}

\alias{RPA.iteration}

\title{Estimating model parameters d and sigma2.}

\description{Finds point estimates of the model parameters d (estimated
true signal underlying probe-level observations), and sigma2
(probe-specific variances).}

\usage{RPA.iteration(S, epsilon = 1e-3, 
                     alpha = NULL, beta = NULL,
                     sigma2.method = "fast", d.method = "fast",
                     maxloop = 1e6)}

\arguments{
  
  \item{S }{Matrix of probe-level observations for a single probeset:
  samples x probes.}
  
  \item{epsilon }{Convergence tolerance. The iteration is deemed
  converged when the change in all parameters is < epsilon.}
  
  \item{alpha, beta }{Priors for inverse Gamma distribution of
  probe-specific variances. Noninformative prior is obtained with alpha,
  beta -> 0.  Not used with sigma2.method 'var'. Scalar alpha and beta
  are specify equal inverse Gamma prior for all probes to regularize the
  solution. The defaults depend on the method.}

  \item{sigma2.method }{

    Optimization method for sigma2 (probe-specific variances).

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

   \item{maxloop }{
		Maximum number of iterations in the estimation process.}
	      
}

\details{Assuming data set S with P observations of signal d with
Gaussian noise that is specific for each observation (specified by a
vector sigma2 of length P), this method gives a point estimate of d
and sigma2. Probe-level variance priors alpha, beta can be used with
sigma2.methods 'robust', 'mode', and 'mean'.  The d.method = "fast" is
the recommended method for point computing point estimates with large
sample size.}

\value{

  A list with the following elements:

  \item{d }{A vector. Estimated 'true' signal underlying the noisy
  probe-level observations.}
  
  \item{sigma2 }{A vector. Estimated variances for each measurement (or
  probe).}
  
}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE 2011.}

\author{Leo Lahti \email{leo.lahti@iki.fi}}

\examples{

## Not run:

## Preprocess probe-level data
## cind determines the 'reference' array
#Smat <- RPA.preprocess(Dilution, cind = 1)

## Pick probe-level data for one probe set
#pmindices <- pmindex(Dilution, "1000_at")[[1]]
#S <- t(Smat$fcmat[pmindices, ])

## RPA with default parameters:
#res <- RPA.iteration(S)

}


\keyword{ methods }
\keyword{ iteration }