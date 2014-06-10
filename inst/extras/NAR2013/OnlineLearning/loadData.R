print("Sample classes")
load("data/classInfo.RData") # classInfo
emat.files <- c(rpa = "data/Lukk-RPA.RData",
	        rpa.online = "~/Rpackages/RPA/OnlineLearning/emat.rpa.online.storage.RData",  
	        rma = "data/Lukk-RMA.RData", 
		frma = "data/Lukk-fRMA.RData", 
		frmab = "data/Lukk-fRMA.batch.RData", 
		mas = "data/Lukk-MAS5.RData")

set.seed(324353)
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
emat <- emat.rpa[, s]
samples <- colnames(emat)

##########################################################

#load("data/Lukk-MAS5-HGU133A_Hs_ENSG.RData")
#load("data/Lukk-fRMA.batch-HGU133A_Hs_ENSG.RData")  
#load("data/Lukk-RMA-HGU133A_Hs_ENSG.RData")
#load("data/Lukk-fRMA-HGU133A_Hs_ENSG.RData")        
#load("data/Lukk-RPA-HGU133A_Hs_ENSG.RData")
#emat.rpa.ENSG
#emat.frma.ENSG
#emat.frma.batch.ENSG
#emat.rma.ENSG
#emat.MAS5.ENSG

