# This file is a part of the RPA program (Robust Probabilistic
# Averaging), see http://www.cis.hut.fi/projects/mi/software/RPA/
#
# Copyright (C) 2008-2010 Leo Lahti (leo.lahti@iki.fi)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License 2 for more details.
# 


initialize.priors <- function (abatch, sets, alpha = 1e-6, beta = 1e-6, d = NULL) {

	# set priors for d and sigma2
	# prior for d omitted in the current package version

	priors <- vector(length = length(sets), mode = "list")

	for (set in sets) {	
		nprobes <- nrow(pm(abatch, set))
		 alphas <- rep.int(alpha, nprobes)
		  betas <- rep.int(beta, nprobes)
	  priors[[set]] <- list(alpha = alphas, beta = betas, d = d) 
	}

	new("rpa.priors", priors)
}

