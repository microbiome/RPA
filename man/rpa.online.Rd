\name{rpa.online}

\Rdversion{1.1}

\alias{rpa.online}

\title{RPA-online for preprocessing very large expression data sets.}

\description{rpa.online is used to preprocess very large expression data
  collections based on a Bayesian hyperparameter update
  procedure. Returns an expressionSet object preprocessed with RPA.}

\usage{rpa.online(cel.path = NULL, cel.files = NULL, sets = NULL, cdf =
                 NULL, bg.method = "rma", priors = list(alpha = 2, beta = 1), epsilon = 1e-2, mc.cores = 1,                 
                 verbose = TRUE, shuffle = TRUE, batch.size = 100,
                 batches = NULL, quantile.basis = NULL, save.batches = TRUE, 
		 save.batches.dir, keep.batch.files = FALSE, unique.run.identifier = NULL)}

\arguments{
  
  \item{cel.path}{Path to CEL file directory}
  
  \item{cel.files}{List of CEL files to preprocess}
  
  \item{sets }{Probesets for which RPA will be computed}

  \item{cdf }{Specify an alternative CDF environment}

  \item{bg.method }{Specify background correction method. See bgcorrect.methods() for options.}
  
  \item{priors }{An 'rpa.priors' object. Can be used to set
		user-specified priors for the model parameters. Not
		used tau2.method = "var". The prior parameters alpha
		and beta are prior parameters for inverse Gamma
		distribution of probe-specific
		variances. Noninformative prior is obtained with
		alpha, beta -> 0.  Not used with tau2.method
		'var'. Scalar alpha and beta specify an identical
		inverse Gamma prior for all probes, which regularizes
		the solution. Can be also specified as lists, each
		element corresponding to one probeset.}

  \item{epsilon }{Convergence tolerance. The iteration is deemed converged 
		  when the change in all parameters is < epsilon.}
  
  \item{mc.cores}{Number of cores for parallel computation}

  \item{verbose}{Print progress information during computation}

  \item{shuffle}{Form random batches}  

  \item{batch.size}{Batch size for online mode (rpa.online); the
complete list of CEL files will be preprocessed in batches with this
size using Bayesian online-updates for probe-specific parameters.}

  \item{batches}{User-defined CEL file batches}

  \item{quantile.basis}{Pre-calculated basis for quantile normalization}

  \item{save.batches}{Save batches to speed up preprocessing.}

  \item{save.batches.dir}{Output directory for temporary batch saves.}

  \item{keep.batch.files}{Logical. Keep (TRUE) or remove (FALSE) the batch files after preprocessing.}

  \item{unique.run.identifier}{Define identifier for this run for naming the temporary batch files. By default, a random id is generated.}

}

\details{RPA preprocessing function. Gives an estimate of the
  probeset-level mean parameter d of the RPA model, and returns these in
  an expressionSet object. The CEL files are handled in batches to
  obtain Bayesian updates for probe-specific hyperpriors; after sweeping
  through the database in batches the results are combined. The online
  mode is useful for preprocessing very large expression data sets where
  ordinary preprocessing algorithms fail, without compromises in
  modelling stage.}

\value{An instance of the 'expressionSet' class.}

\references{See citation("RPA").}

\author{Leo Lahti \email{leo.lahti@iki.fi}}

\note{Development version}

\seealso{rpa, AffyBatch, ExpressionSet}

\examples{

# Not run:
# eset <- rpa.online(cel.file.path)

}

\keyword{ methods }

