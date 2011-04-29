setClass("rpa", 	 
	 representation(d = "array", mu.real = "list", sigma2 = "list", affinity = "list", cind = "numeric", sets = "character", data = "array", cdf = "character", abatch = "AffyBatch"), 
	 contains = "list")

setClass("rpa.fit", 	 
	 representation(mu = "array", mu.real = "numeric", sigma2 = "list", affinity = "list", data = "array"), 
	 contains = "list")

setClass("rpa.list",
         representation(d = "array", mu.real = "numeric", sigma2 = "list", affinity = "list", cind = "numeric", set = "character", data="array"), 
	 contains = "list")

setClass("rpa.priors",
	 representation(alpha = "list", beta = "list", d = "list"),
	 contains = "list")
