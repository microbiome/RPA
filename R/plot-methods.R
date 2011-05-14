#setMethod(f = "plot", signature("rpa.list"),
#        function(x, y){
#        }
#)


plot.rpa.list <- function (x, y, ...) {

      set <- x$set # probe set
       mu <- x$d # 'Absolute' expression level
       sd <- sqrt(x$sigma2) # probe standard deviation	 

      par(mfrow = c(2,1))
      barplot(sd, main = paste(set,"/ Probe noise (standard deviation)"),
              ylab = "Standard deviation", xlab = "Probes", names.arg = 1:length(sd), las = 2)

      barplot(mu, main = paste(set,"/ expression "),
              ylab = "Expression signal", xlab = "Arrays", las = 1) 

}

plot.rpa <- function (x, y = NULL, set, highlight.probes = NULL, pcol = "darkgrey", mucol = "black", ecol = "red", cex.lab = 1.5, cex.axis = 1, cex.names = 1, cex.main = 1, external.signal = NULL, main = NULL, plots = "all", ...) {
                                      
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
    main <- paste("Probe signals and the summary estimate (", set, ")", sep = "")
  }
  rpa.plot(dat, rpa.fit.object, cex.lab = cex.lab, cex.axis = cex.axis, cex.names = cex.names, cex.main = cex.main, main = main, plots = plots, external.signal = external.signal, pcol = pcol, ecol = ecol, mucol = mucol)

  dat
}




