
#                                        
# This file is a part of the RPA program (Robust Probabilistic
# Averaging), see http://www.cis.hut.fi/projects/mi/software/RPA/
#
# Copyright (C) 2008-2010 Leo Lahti (leo.lahti@iki.fi)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License 2 for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#
  
rpa <- function (abatch,
                    sets = NULL,
                  myseed = 101, 
                  priors = NULL,
                 epsilon = 1e-2, 
                    cind = 1,
           sigma2.method = "robust",
                d.method = "fast",
                 verbose = FALSE,
               bg.method = "rma",
    normalization.method = "quantiles.robust",
                     cdf = NULL,
                   alpha = NULL,
                    beta = NULL,
                 exclude.reference.array = FALSE) 
{

#
# RPA for preprocessing only
#

  #Set random seed
  set.seed( myseed )

  # Set alternative CDF environment if given
  if (!is.null(cdf)) { abatch@cdfName <- cdf }

  # Preprocessing
  preproc <- RPA.preprocess(abatch, cind, bg.method, normalization.method, cdf)	
  
  # Pick the necessary objects
      fcmat <- preproc$fcmat      # probe-level differential expressions
       cind <- preproc$cind       # index of the control array
   set.inds <- preproc$set.inds   

  #################################################################

  # ESTIMATE PROBE RELIABILITY AND DIFFERENTIAL GENE EXPRESSION

  #################################################################

  message("Calculating Expression")
  
  # Number of arrays except control
  T <- ncol( fcmat )

  # Check names and number for the investigated probesets
  # if my.sets not specified, take all sets in abatch
  if (is.null(sets)) { sets <- geneNames(abatch) } 
  Nsets <- length( sets ) 

  ## Matrices to store results
  d.results <- array(NA, dim = c(Nsets, T))
  rownames(d.results) <- sets
  colnames(d.results) <- colnames(fcmat)

  if (!is.null(priors) && (!is.null(alpha) || !is.null(beta))) {
    stop("Specify either priors OR alpha, beta- both cannot be specified at the same time!")
  }

  for (i in 1:Nsets) {

    set <- sets[[i]]
  
    if (verbose) {message(paste("Computing probeset", set, ":", i, "/", Nsets, "...\n"))}

    # Find probe (pm) indices for this set
    pmindices <- set.inds[[set]] #pmindex(abatch, set)[[1]]
  
    # Number of probes for this probeset
    P <- length(pmindices)

    # Get chips x probes matrix of probe-wise fold-changes
    #S <- t(fcmat[pmindices, ])
    S <- t(fcmat[pmindices, ])
    
    # Pick the priors for this set
    if (!is.null(priors)) {
      alpha <- priors[[set]]$alpha 
      beta  <- priors[[set]]$beta
    } else {}
        
    # Calculate RPA 
    res <- RPA.iteration(S, epsilon, alpha, beta, sigma2.method, d.method)
                                                  
    # Store the results (only mean parameter)
    d.results[i, ] <- res$d
    
  }

  # Excluding reference array (if set so). All values are zero since #
  #the method gives differential expression of the reference array #
  #against itself. Shifting the reference this way does not affect the
  #relative differences # between the arrays.

  if (!exclude.reference.array) {
    d.ref <- rep.int(0, nrow(d.results))
    dat <- cbind(d.ref, d.results)
    colnames(dat) <- c(sampleNames(abatch@phenoData)[[cind]], colnames(d.results))
    rownames(dat) <- rownames(d.results)

    # Order samples to the same order as in original affybatch
    dat <- dat[,sampleNames(abatch@phenoData)]
    
  } else {dat <- d.results}
  
  
  # Coerce expression values in the rpa object into an ExpressionSet object
  # and return expression set
  return( new("ExpressionSet", assayData = list(exprs = dat)) )
  
}

