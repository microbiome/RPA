 #.First.lib <- function (libname, pkgname, where) {
	   #if (!require(methods)) stop("We require methods for package RPA")
	   #where <- match(paste("package:", pkgname, sep=""), search())
	   #.initRPAmethods(where)
	   #.initRPAmethods()
#	   NULL
#}






setClass("rpa", 	 
	 #representation(d="matrix",sigma2="list",cind="numeric",sets="character"), 
	 contains = "list")

setClass("rpa.priors",
	 #representation(alpha="list",beta="list"),
	 contains = "list")

setClass("rpa.list", 
	 #representation(d="matrix",sigma2="list",cind="numeric",set="character"), 
	 contains = "list")



#setGeneric("[",function(x,i,j,drop){standardGeneric("[")})
setMethod("[", "rpa",
   function(x, i, j, ... , drop=TRUE) {
      if (typeof(i) == "character") {i=which(x$sets == i)}
      #new("rpa.list", list(d = x@d[i,], sigma2=x@sigma2[[i]], cind = x@cind, set = x@sets[[i]]))
      new("rpa.list", list(d = x$d[i,], sigma2=x$sigma2[[i]], cind = x$cind, set = x$sets[[i]]))
    }
)
			   
			   
#setGeneric("[[",function(x,i,j,drop){standardGeneric("[[")})
setMethod("[[", "rpa",
   function(x, i, j="missing", ..., exact=TRUE) {
      if (typeof(i) == "character"){i=which(x$sets == i)}
      new("rpa.list", list(d = x$d[i,], sigma2=x$sigma2[[i]], cind = x$cind, set = x$sets[[i]]))
   }
)


#setGeneric("show",function(object){standardGeneric("show")})
setMethod(show, "rpa",function(object) {cat("rpa object\n")}) #{show.rpa()})
#show.rpa <- function (...) {cat("rpa object \n")}

plot.rpa.list <- function (x, y, ...) {

      set <- x$set
      d <- x$d
      s2 <- x$sigma2		  

      par(mfrow=c(2,1))
      barplot(s2,main=paste(set,"/ Probe-specific variances (sigma2)"),ylab="Variance",xlab="Probes",
              names.arg=paste(set,seq(length(s2))),las=2)
      barplot(d,main=paste(set,"/ Differential expression estimate (d)"),
              ylab="Dif. exp. signal", xlab="Arrays",las=2)	      

}









