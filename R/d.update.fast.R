# This file is a part of the RPA program
# (Robust Probabilistic Averaging) 
# http://bioconductor.org/packages/release/bioc/html/RPA.html

# Copyright (C) 2008-2011 Leo Lahti <leo.lahti@iki.fi>. All rights reserved.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the FreeBSD License.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


d.update.fast <- function (St, s2) {

  # St: probes x arrays  

  # With large sample sizes when T -> Inf
  # the d converges to the weighted mean 
  # over the probes, weighted by probe variances	 
  colSums(St / s2) / sum(1 / s2)
}


#require(compiler)
d.update.fast.c <- cmpfun(d.update.fast) 
