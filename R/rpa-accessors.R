setMethod(f = "[", signature("rpa"), 
   definition = (function(x, i, j, ... , drop=TRUE) {
      if (typeof(i) == "character") { i = which(x$sets == i) }
      new("rpa.list", list(d = x$d[i,], sigma2 = x$sigma2[[i]], cind = x$cind, set = x$sets[[i]]),)
    })
)

setMethod(f = "[[", signature("rpa"),
   definition = (function(x, i, j = "missing", ..., exact = TRUE) {
      if (typeof(i) == "character"){i = which(x$sets == i)}
      new("rpa.list", list(d = x$d[i,], sigma2 = x$sigma2[[i]], cind = x$cind, set = x$sets[[i]]))
   })
)

#setReplaceMethod(f = "[[", signature("rpa"),
#   definition = (function(x, i, j, value) {
#      x@models[[i]] <- value                                        
#      return(x)                                        
#      }                               
#))

