# 26.5.2012 Summaries / observations for the manuscript

# Final preprocessings using RPA-online-storage mode
source("Lukk-RPA-ENSG-20120528.R")
source("Lukk-RPA-20120528.R")

# hyperparameter distributions & evolution
source("investigate.hyperparameters.R")

# Run time test -> OK 6.6.2012: see runtimes.eps
source("runtimetest.R") # effect of mc.cores & batch.size?
# -> running time plot for collection of 20-30k CEL files (with multiple cores?); mention speedups are possible

# TODO Experiments
# - RPA utilizes information from the complete data set; show
evolution of probe effects (standard deviation) and claim that there
are many probes that cannot be deemed converged even after; 5000
samples; how is fRMA with 1000 samples -> suggest a better minimum sample size for robust estimation : 3000?
# - include probe effect data tables or RData files as supplementary information (affinities & variances)
# - Note: estimates of d would converge more rapidly than probe-specific variances: show results?
# - Note: probe variances shown to match with known error sources already in TCBB/IEEE. Now we also have affinities. Would be cool the check these in more detail but out of scope,

# TODO Package devel
#- report probe effect tables also in the rpa and rpa.online outputs?
#- For speeupds, it is possible to use RPA in the same way as fRMA - estimate probe effects only from a subset: provide also this option in rpa implementation?

# RMA, RPA, RPA-Online: comparisons
# source("online.version.R") # -> variant.correlations.tab
# This confirms that RPA-online with batch storage produces essentially the same results as the standard RPA 
# (correlation > 0.996 for 300 random CEL files and batch of 50)
# For some reason, the correlation to RMA is also best with the RPA-online storage variant (r = 0.965)
# TODO: marginalizing out d in the probe effect estimation by sampling d for each sample, and 
#    solving taus based on this, then using the mode of tau as a robust estimate of the posterior

############################################

# OTHER USEFUL STUFF MISC

# Test hyperparameter estimation with toydata. 
source("test2.R")
source("example2.R")
source("known.abatch.R") # inserting toy probeset in real data

# Run the pipeline in steps
source("EBI-Lukk-AltCDF-RPA.R")
source("test.ebi.R")

# Calculates correlation of the estimated
# signal and real signal in online-learning and alternatives.
# When batch size grows in online-learning ("step")
# the result converges to RPA calculated from the complete data
# When batch size is small compared to the data, 
# online learning is not better than median
source("online-tests.R")
source("online.test.R")
source("online.test2.R") # OUTOA: batch-filujen kanssa saa paremmat korrelaatiot RPA vs RPA-online? 

# Compare quantile-normalized data from the batches to quantiles in standard preprocessing
# Assumes the batches have been calculated elsewhere
source("validate.normalization.R")

# testing with SNP probes
source("altcdf.R")
source("altcdf-ENSEMBL.R")

# parallelization from Triton
source("parallel.R")
