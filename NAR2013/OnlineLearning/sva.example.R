
     ## Load data
       library(bladderbatch)
       data(bladderdata)
       
       ## Obtain phenotypic data
       pheno = pData(bladderEset)
       edata = exprs(bladderEset)
       batch = pheno$batch
       mod = model.matrix(~as.factor(cancer), data=pheno)
       mod0 = model.matrix(~1, data=pheno)
       
       ## Construct the surrogate variables 
       svaobj <- sva(edata,mod,mod0,method="irw")
      
       ## Include them in a downstream analysis
       mod.sv <- cbind(mod,svaobj$sv)
       mod0.sv <- cbind(mod0,svaobj$sv)
       adjusted.pvals <- f.pvalue(dat,mod.sv,mod0.sv)
