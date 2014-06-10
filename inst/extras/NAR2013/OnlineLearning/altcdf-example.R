library(hgu133plus2hsensgcdf) # biocLite("hgu133plus2hsensgcdf")
library(CustomCDF) # http://arrayanalysis.mbni.med.umich.edu/MBNIUM.html
cels <- list.celfiles("/my.cel.path/", full.names = T)
abatch <- ReadAffy(filenames = cels)
abatch@cdfName <- "HGU133Plus2_Hs_ENSG"                                         
