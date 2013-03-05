# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2012 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


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

