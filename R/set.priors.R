# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


set.priors <- function (abatch, set, 
			alpha, beta, priors = NULL, 
			alpha.template = 1e-6, beta.template = 1e-6) {

  if (is.null(priors)) {
    priors <- initialize.priors(abatch,
                                 sets = geneNames(abatch),
                                alpha = alpha.template,
                                 beta = beta.template)
  }

  priors[[set]]$alpha <- alpha
  priors[[set]]$beta <- beta

  priors
}
