\name{RPA.iteration}

\Rdversion{1.1}

\alias{RPA.iteration}

\title{Finding a mode for the model parameters d and sigma2.}

\description{Finds a mode for the model parameters d (estimated true
signal underlying probe-level observations), and sigma2 (probe-specific
variances).}

\usage{RPA.iteration(S, sigma2.guess, epsilon = 10^(-3), 
                     alpha = NULL, beta = NULL,
                     sigma2.method = "var", d.method = "fast")}

\arguments{
  
  \item{S}{Matrix of probe-level observations for a single probeset:
  samples x probes.}
  
  \item{sigma2.guess}{Initial values for probe-specific variances.}
  
  \item{epsilon}{Convergence tolerance. The iteration is deemed
  converged when the change in all parameters is < epsilon.}
  
  \item{alpha, beta}{Vectors that specify priors for inverse Gamma
  distribution of probe-specific variances. Noninformative prior is
  obtained with alpha, beta -> 0.  Used only when sigma2.method =
  "basic".}

  \item{sigma2.method}{
    Optimization method for sigma2 (probe-specific variances).
	
        "basic": optimization using user-specified alpha, beta priors.
		 Default: alpha, beta = 1e-6.
		 
        "var": utilizes the fact that the cost function converges to
	       variance with large sample sizes. Default method.
	
      }

  \item{d.method}{
    Method to optimize d.

        "basic": optimization scheme to find a mode used in Lahti et
        	 al. TCBB/IEEE; relatively slow; this is the preferred 
		 method with small sample sizes.

        "fast": weighted mean over the probes, weighted by probe
                variances The solution converges to this with large
                sample size. Default method.

      }
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

  A list with the following elements:

  \item{d }{A vector. Estimated 'true' signal underlying the noisy
  probe-level observations.}
  
  \item{sigma2 }{A vector. Estimated variances for each measurement (or
  probe).}
  
}

\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE, to appear. See
http://www.cis.hut.fi/projects/mi/software/RPA/}

\author{Leo Lahti <leo.lahti@tkk.fi>}

\examples{

require(affydata)
data(Dilution) # load affybatch

# select array against which differential expression is computed the
# method has appeared to be robust for the choice of cind
cind <- 1

# Preprocess probe-level data
Smat <- RPA.preprocess(Dilution,cind)

# Pick probe-level data for one probe set
set <- "1000_at"
pmindices <- pmindex(Dilution, set)[[1]]
S <- t(Smat$fcmat[pmindices, ])

# Set initial values for probe-specific variances
sigma2.guess <- rep(10,ncol(S))

# assume non-informative prior and large sample size
# and se fast optimization
res <- RPA.iteration(S, sigma2.guess,
                     epsilon = 10^(-3),
                     sigma2.method = "var",
                     d.method = "fast")

par(mfrow=c(1,2))

barplot(res$d,names.arg=paste("Array",seq(ncol(exprs(Dilution)))[-cind]),
        main="d",xlab="Arrays/samples",ylab="Differential expression",las=2)

barplot(res$sigma2,las=2,names.arg=paste(set,seq(ncol(S))),
        main="sigma2",xlab="Probes",ylab="Variance",las=2)

# Generate toy data where the ground truth is known
# 300 arrays x 11 probes

set.seed(24)
d.true <- rnorm(300,0,4)
sigma2.true <- c(1,1,1,1,2,2,3,3,3,4,6)
S <- NULL

for (s2 in sigma2.true) {
   S <- cbind(S,rnorm(length(d.true),mean=d.true,sd=sqrt(s2)))
}

# Initial guess for sigma2
sigma2.guess <- rep(1,ncol(S))

# assuming non-informative priors and large sample size
res <- RPA.iteration(S, sigma2.guess,
                     epsilon = 10^(-3),
                     sigma2.method = "var",
                     d.method = "fast")

barplot(res$sigma2,ylab="Estimated var", xlab="Probes",main="sigma2")}

\keyword{ iteration }
\keyword{ methods }
