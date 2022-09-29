print("Sample classes")
load("data/classInfo.RData") # classInfo
emat.files <- c(rpa = "data/Lukk-RPA.RData",
	        rpa.online = "~/Rpackages/RPA/OnlineLearning/emat.rpa.online.RData",  
	        rma = "data/Lukk-RMA.RData", 
		frma = "data/Lukk-fRMA.RData", 
		frmab = "data/Lukk-fRMA.batch.RData", 
		mas = "data/Lukk-MAS5.RData")

set.seed(324353)
#rs <- sample(names(classInfo), 100) # random subset of the data

rs <- names(classInfo)
classInfo <- classInfo[rs]

# Remove smallest classes
tab <- table(classInfo)
class.minsize <- 0
selected.classes <- names(which(tab >= class.minsize) )
s <- names(classInfo)[which(classInfo %in% selected.classes)]
cl <- classInfo[s]

# Pick one of the data matrices
method <- "rpa" # names(emat.files)[[1]]
load(emat.files[[method]])
if (method == "rpa") { emat <- emat.rpa[, s]}
if (method == "rma") { emat <- emat.rma[, s]}
