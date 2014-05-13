\name{updating.hyperparameters}
\alias{updating.hyperparameters}
\title{updating.hyperparameters
Hyperparameter update}
\usage{
  updating.hyperparameters(q, set.inds, verbose,
    mc.cores = 1, alpha, betas, epsilon)
}
\arguments{
  \item{q}{probes x samples matrix}

  \item{set.inds}{Probe set indices}

  \item{verbose}{Print progress information}

  \item{mc.cores}{Number of cores for parallel computation}

  \item{alpha}{alpha hyperparameter}

  \item{betas}{beta hyperparameters}

  \item{epsilon}{Convergence parameter}
}
\value{
  List with the following elements: alpha, betas, s2s
  (variances)
}
\description{
  updating.hyperparameters Hyperparameter update
}
\examples{
#
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("RPA")
}
\keyword{utilities}

