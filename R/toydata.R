sample.probeset <- function (P = 10, n = 20, shape = 1, scale = 1, mu.real = 2) {

  # Generating toy data
  # P: number of probes (observations)
  # n: number of samples
  # shape, scale: parameters of the (inverse) Gamma conjugate prior
 
  # Sample probe-specific variances from 
  # inverse Gamma distribution a.k.a scaled inverse chi-square
  # (this is conjugate prior for variances, 
  probe.variance <- 1/rgamma(P, shape = shape, scale = scale) # invgam ~ 1/gam
  
  # RPA model assumes that affinities are also distributed as 
  # N(0, sigmaj2)
  probe.affinity <- rnorm(P, sd = sqrt(probe.variance))

  # Generate real signal shape (d)
  d <- rnorm(n)

  dat <- array(NA, dim = c(P, n))
  for (p in 1:P) {
    dat[p, ] <- rnorm(n, 
    	              mean = d + mu.real + probe.affinity[[p]], 
                        sd = sqrt(probe.variance[[p]]))
  }

  list(dat = dat, variance = probe.variance, affinity = probe.affinity, d = d, mu.real = mu.real)

}

