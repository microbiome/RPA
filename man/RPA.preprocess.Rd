\name{RPA.preprocess}

\Rdversion{1.1}

\alias{RPA.preprocess}

\title{Preprocess AffyBatch object for RPA.}

\description{Quantile normalizes and log2-transforms the probe-level raw
data in abatch. Then probe-level differential expression is computed
between the specified 'control' array (cind) and the other arrays.}

\usage{RPA.preprocess(abatch, cind)}

\arguments{

  \item{abatch}{An AffyBatch object.}
  
  \item{cind}{Specifies which of the arrays in abatch is used as a
  control in computing probe-level differential expression.}
  
}

\details{The probe-specific variance estimates are robust against the
choice of cind; cind can be selected at random.}

\value{
  
  \item{fcmat }{Probes x arrays preprocessed differential expression matrix.}

  \item{cind }{Specifies which array in abatch was selected as a control
  for computing probe-level differential expression.}
  
}


\references{Probabilistic Analysis of Probe Reliability in Differential
Gene Expression Studies with Short Oligonucleotide Arrays.  Lahti et
al., TCBB/IEEE, to appear. See
http://www.cis.hut.fi/projects/mi/software/RPA/ }

\author{Leo Lahti <leo.lahti@tkk.fi>}

\seealso{AffyBatch}

\keyword{ methods }

