library(CustomCDF)
require(RPA)

set.seed(11122) 
cels <- list.celfiles("/data/GSE3526/CEL/", full.names = T)[1:3]

abatch <- ReadAffy(filenames = cels)
set <- "1494_f_at"
print(paste("Number of probes for probeset ", set, ":", nrow(pm(abatch, set))))
print(pmindex(abatch, set))

xy <- getsnpprobe(chip = "HG-U133_Plus_2")
report <- removeprobe(abatch, xyMatrix = xy)
print(paste("Number of probes for probeset ", set, ":", nrow(pm(abatch, set))))
print(dim(xy))
print(pmindex(abatch, set))

abatch@cdfName <- "HGU133Plus2_Hs_ENSG"   

#removeprobe(object, xyMatrix=NULL

#########################################
#eset <- rma(abatch)
#set <- report$removedprobe$probeset[[1]]
#barplot(rbind(exprs(eset.filtered)[set,], exprs(eset)[set,]), beside = TRUE)

#xy <- getsnpprobe(chip = abatch@cdfName)
#xy <- getsnpprobe(chip = "HGU133Plus2")
#xy<- getsnpprobe(chip='HGU133Plus2', criteria='het > 0.5')
#cdf <- "HGU133Plus2_Hs_ENSG"                                        
#if (!is.null(cdf)) { abatch.orig@cdfName <- cdf }

# loc > 0 and loc < 10
# het > 0
# rsid in (1,3,4,5)
# loc   Location of allele on the probe, upper value is probe's length,
#       lower limit could be negative if the allele is longer than probe
# het	Heterozygosity value
# rsid	Snp id



    
