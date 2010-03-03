setClass("rpa", 	 
	 representation(d = "array", sigma2 = "list", cind = "numeric", sets = "character", data = "array", cdf = "character", abatch = "AffyBatch"), 
	 contains = "list")

setClass("rpa.list",
         representation(d = "array", sigma2 = "list", cind = "numeric", set = "character", data="array"), 
	 contains = "list")

setClass("rpa.priors",
	 representation(alpha = "list", beta = "list", d = "list"),
	 contains = "list")
