RPA.dcost <- function (d, tau2, S) {

  #
  # S : observed probe-level signals, T x P i.e. arrays x probes
  #
  # d : assumed "real" signal for which this is a cost function
  #
  # tau2 : probe-specific variances

  R <- S - d

  -sum((1/(2*tau2))*((colSums(R)^2)/(nrow(R) + 1) - colSums(R^2)))

}

