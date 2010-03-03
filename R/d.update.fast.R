#
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

d.update.fast <- function (St, sigma2) {
  # With large sample sizes when T -> Inf
  # the d converges to the weighted mean 
  # over the probes, weighted by probe variances	 
  colSums(St / sigma2) / sum(1 / sigma2)
}
