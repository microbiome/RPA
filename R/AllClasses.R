setClass("rpa.list",
         representation(d = "array", mu.real = "numeric", tau2 = "list", affinity = "list", cind = "numeric", set = "character", data="array"), 
	 contains = "list")

