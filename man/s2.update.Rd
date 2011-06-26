\name{s2.update}
\alias{s2.update}
\title{Variance update.}
\description{Update variance parameters in the Bayesian RPA model with an EM-type optimization procedure and potentially informative priors.}
\usage{
s2.update(dat, alpha = 1e-2, beta = 1e-2, s2.init = NULL, th = 0.01)
}

\arguments{
  \item{dat}{A probes x samples matrix (probeset).}
  \item{alpha, beta}{Shape and scale parameters of the inverse gamma prior for the variances.}
  \item{s2.init}{Initial values for the variances for optimization.}
  \item{th}{Optimization convergence threshold.}
}

\details{Updates the variance parameters in the Bayesian RPA model
with a (generalized) EM-type optimization procedure and potentially
informative priors. The variances are updated by their mode at each
step until convergence. Priors (alpha, beta) are taken into account if
provided.}

\value{A vector of variance parameters, one for each probe.}

\references{See citation("RPA").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
%\note{}

\seealso{rpa.online}
\examples{
# dat: probes x samples matrix (probeset)
# s2 <- s2.update(dat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ iteration }


