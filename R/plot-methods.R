#setMethod(f = "plot", signature("rpa.list"),
#        function(x, y){
#        }
#)


plot.rpa.list <- function (x, y, ...) {

      set <- x$set
        d <- x$mu
       s2 <- x$sigma2		 

      par(mfrow = c(2,1))
      barplot(s2, main = paste(set,"/ Probe-specific variances (sigma2)"),
              ylab = "Variance", xlab = "Probes", names.arg = 1:length(s2), las = 2)              
      #barplot(d, main = paste(set,"/ Differential expression estimate (d)"),
      barplot(d, main = paste(set,"/ Expression estimate"),
              ylab = "Signal", xlab = "Arrays", las = 2)	      

}

plot.rpa <- function (x, y = NULL, set, highlight.probes = NULL, pcol = "darkgrey", mucol = "black", ecol = "red", cex.lab = 1.5, cex.axis = 1, external.signal = NULL, main = NULL, ...) {
                                      
  # x: 'rpa' object to plot  
  # set: probeset

  # get the associated affybatch
  abatch <- x$abatch
  
  # Use alternative CDF environment if given                                        
  if (!is.null(x$cdf)) { abatch@cdfName <- x$cdf }

  # Find probe (pm) indices for this set
  pmindices <- which(probeNames(abatch) == set)
  
     d <- x[[set]]$d
    s2 <- x[[set]]$sigma2
  cind <- x[[set]]$cind
   dat <- x$data[pmindices,]
  mu.real <- x[[set]]$mu.real
  affinity <- x[[set]]$affinity
  
  rpa.fit.object <- new("rpa.fit", list(mu = d, mu.real = mu.real, sigma2 = s2, affinity = affinity, data = dat))

  if (is.null(main)) {
    main <- paste(set, "/ Probe-level signals and the summary estimate")
  }
  rpa.plot(dat, rpa.fit.object, cex.lab = cex.lab, cex.axis = cex.axis, main = main, plots = "data", external.signal = external.signal, pcol = pcol, ecol = ecol, mucol = mucol)

  dat
}




