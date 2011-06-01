\name{RPA.preprocess}

\Rdversion{1.1}

\alias{RPA.preprocess}

\title{Preprocess AffyBatch object for RPA.}

\description{Background correction, quantile normalization and
  log2-transformation for probe-level raw data in abatch. Then
  probe-level differential expression is computed between the specified
  'reference' array (cind) and the other arrays.}

\usage{RPA.preprocess(abatch, cind = 1, bg.method = "rma",
                      normalization.method = "quantiles.robust", cdf = NULL, quantile.n = 50)}
                            

\arguments{

  \item{abatch }{An AffyBatch object.}
  
  \item{cind }{Specify which of the arrays in abatch is used as a
  reference for computing probe-level differential expression.}

  \item{bg.method }{Specify background correction method. See
  bgcorrect.methods(abatch) for options.}
  
  \item{normalization.method }{Specify normalization method. See
  normalize.methods(abatch) for options. For memory-efficient online
  version, use "quantiles.online".}

  \item{cdf }{The CDF environment used in the analysis.}
  
  \item{quantile.n }{Number of samples to be used in quantiles.online normalization.}
}

\details{Probe-specific variance estimates are robust against the choice
of reference array.}

\value{
  
  \item{fcmat }{Probes x arrays preprocessed differential expression matrix.}

  \item{cind }{Specifies which array in abatch was selected as a reference
               in calculating probe-level differential expression.}

  \item{cdf }{The CDF environment used in the analysis.}

  \item{set.inds }{Indices for probes in each probeset, corresponding to
                   the rows of fcmat.}

}

\references{See citation("RPA").}

\author{Leo Lahti <leo.lahti@iki.fi>}

\seealso{AffyBatch}

\keyword{ methods }

