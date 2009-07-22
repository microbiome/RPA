\name{RPA.pointestimate}

\Rdversion{1.1}

\alias{RPA.pointestimate}

\title{Computing point estimate for the model parameters for all probe
sets.}

\description{Computes point estimate}

\usage{RPA.pointestimate(abatch, sets = NULL, myseed = 101, 
			 priors = NULL, 
			 epsilon = 10^(-2), cind = FALSE, 
			 sigma2.method = "var",
			 d.method = "fast",
			 verbose = TRUE)}

\arguments{

  \item{abatch }{An AffyBatch object.}

  \item{sets }{Specifies the probesets for which RPA estimates will be
  computed. Default: all probe sets.}

  \item{myseed }{Specifies the random seed.}

  \item{priors }{An 'rpa.priors' object. Can be used to set
		user-specified priors for the model parameters. Used
		only with sigma2.method = "basic".}
  
  \item{epsilon }{Convergence tolerance. The iteration is deemed
  converged when the change in all parameters is < epsilon.}
  
  \item{cind }{Specifies which of the arrays in abatch is used as a
  control in computing probe-level differential expression.}


  \item{sigma2.method }{
  	     
	   Optimization method for sigma2 (probe-specific variances).
   
              "basic": optimization using user-specified alpha, beta priors.
                       Default: alpha, beta = 1e-6.

	      "var": utilizes the fact that the cost function converges to
                     variance with large sample sizes. Default method. 
	      }
		     
   \item{d.method }{
	   Method to optimize d.
			   
	      "basic": Finds a mode by directly optimizing the cost
	      	       function of the probabilistic model. 
				  
	      "fast": weighted mean over the probes, weighted by
		      inverse probe-specific variances. The solution
		      converges to this with large sample size.
		      Default method.}
	      
   \item{verbose }{Print progress information during computation. Default: TRUE.}

}

\details{Assuming data set S with P observations of signal d with
Gaussian noise that is specific for each observation (specified by a
vector sigma2 of length P), this method gives a point estimate of d
and sigma2. Note that probe-level variance priors alpha, beta are used
only when sigma2.method = "basic". The sigma2.method = "var" assumes
non-informative priors. The d.method = "fast" is the preferred method
for point computing point estimates when sample size is large. It
computes the average over probe-level observations, weighted by the
inverse probe-specific variances, and is expected to be more robust
and faster than d.method="basic" that finds point estimate for d by
directly optimizing the posterior distribution.}

\value{

  An instance of the 'rpa' class. This is an extended list containing
  the following elements:

  \item{d }{A matrix of probesets x arrays. Specifies the estimated
  'true' underlying differential gene expression signal over the arrays
  (vs. the control array 'cind') for each investigated probeset. Note
  that the control array is not included.}
  
  \item{sigma2 }{A list. Each element corresponds to a probeset, and
  contains a vector that gives the estimated variance for each probe in
  that probeset.}
  
  \item{cind }{Specifies which of the arrays in the abatch (the
	     affybatch object to be analyzed) was used as a control in
	     computing probe-level differential expression.}

  \item{sets }{A character vector listing the investigated probesets.}
}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE, to appear. See
http://www.cis.hut.fi/projects/mi/software/RPA/ }

\author{Leo Lahti <leo.lahti@tkk.fi>}

\note{With large sample size, sigma2.method="var" and d.method="fast"
are recommended. With small sample size and informative prior,
sigma2.method="basic" and d.method="basic" may be preferable.}

\seealso{RPA.iteration, initialize.priors, AffyBatch}

\examples{

# Load example data set
require(affydata)
data(Dilution)

## Run RPA analysis	
## Compute RPA for the whole data set	
## Slow, not executed here	
##rpa.results <- RPA.pointestimate(Dilution)
	
# Compute RPA for specific probesets only	
sets <- geneNames(Dilution)[1:5]	
rpa.results <- RPA.pointestimate(Dilution,sets)

# Visualize the results for one of the probe sets
plot(rpa.results[sets[[1]],])
}

\keyword{ methods }

