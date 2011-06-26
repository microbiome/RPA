\name{hyperparameter.update}
\alias{hyperparameter.update}
\title{Update hyperparameters.}
\description{Update shape (alpha) and scale (beta) parameters of the inverse gamma distribution.}
\usage{
hyperparameter.update(dat, alpha, beta, th = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{A probes x samples matrix (probeset).}
  \item{alpha, beta}{Shape and scale parameters of inverse gamma density for the probe variances.}
  \item{th}{Convergence threshold.}
}
\details{Shape update: alpha <- alpha + T/2; Scale update: beta <- alpha * s2 where s2 is the updated variance for each probe (the mode of variances is given by beta/alpha). The variances (s2) are updated by EM type algorithm, see s2.update.}
\value{ A list with elements alpha, beta (corresponding to the shape and scale parameters of inverse gamma distribution, respectively).}

\references{See citation("RPA").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}

\seealso{s2.update, rpa.online}
\examples{

# Generate and fit toydata, learn hyperparameters
set.seed(11122)
P <- 11   # number of probes
N <- 5000 # number of arrays
real <- sample.probeset(P = P, n = N, shape = 3, scale = 1, mu.real = 4)
dat <- real$dat # probes x samples

# Set priors
alpha <- 1e-2
beta  <- rep(1e-2, P)

# Operate in batches
step <- 1000
for (ni in seq(1, N, step)) {
  batch <- ni:(ni+step-1)  
  hp <- hyperparameter.update(dat[,batch], alpha, beta, th = 1e-2)
  alpha <- hp$alpha
  beta <- hp$beta
}

# Final variance estimate
s2 <- beta/alpha

# Compare real and estimated variances
plot(sqrt(real$variance), sqrt(s2), main = cor(sqrt(real$variance), sqrt(s2))); abline(0,1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ iteration }

