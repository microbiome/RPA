library(CustomCDF)
require(RPA)

# Load altCDF for ENSG and ENST
#if (!requireNamespace("BiocManager", quietly=TRUE))
    #install.packages("BiocManager")
#BiocManager::install("hgu133ahsensgcdf")
#BiocManager::install("hgu133ahsenstcdf")

# Read Affybatch
set.seed(11122) 
cels <- list.celfiles("/share/mi/data/GSE3307/", full.names = T)[1:3]
abatch <- ReadAffy(filenames = cels)

# Optional: Change CDF
#abatch@cdfName <- "HGU133A_Hs_ENSG"
#abatch@cdfName <- "HGU133A_Hs_ENST"                                         

# Total number of probes
np1 <- length(unlist(pmindex(abatch)))
print(paste("Total probes:", np1))

# Check probes with known SNPs
xy <- getsnpprobe(chip = "HG-U133A")

# Remove SNP probes from affybatch
report <- removeprobe(abatch, xyMatrix = xy)

# Check that SNP probes were really removed 
# (abatch now contains less probes than above)
np2 <- length(unlist(pmindex(abatch)))
print(paste("Total probes:", np2))

# How many SNP probes there were?
print(paste("Difference: ", np1 - np2))
print( c(np1, np2, np1-np2))

# With original affy sets (cdfName not set)
#> 247965 - 245657
#[1] 2308
# length(featureNames(abatch))
# ~ 22 283

# ENSG
#> c(np1, np2, np1-np2)
#[1] 171666 171651     15
# length(featureNames(abatch))
#[1] 11869

# ENST
#> print( c(np1, np2, np1-np2))
#[1] 701051 696521   4530
#length(featureNames(abatch))
# > 51956

#> dim(xy)
#[1] 2681    2

In total 2681 SNPs

     Probesets SNP probes 
AFFY 22,283    2308
ENSG 11,869    15
ENST 51,956    4530
