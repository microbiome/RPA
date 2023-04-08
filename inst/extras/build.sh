#~/bin/R-patched/bin/R CMD BATCH document.R
~/bin/R-4.2.2/bin/R CMD build ../../ --resave-data #--no-examples  --no-build-vignettes 
~/bin/R-4.2.2/bin/R CMD check RPA_1.55.2.tar.gz #--no-build-vignettes --no-examples
#~/bin/R-4.2.2/bin/R CMD check --as-cran RPA_1.55.2.tar.gz
~/bin/R-4.2.2/bin/R CMD BiocCheck RPA_1.55.2.tar.gz
~/bin/R-4.2.2/bin/R CMD INSTALL RPA_1.55.2.tar.gz 
#~/bin/R-4.1.0/bin/R CMD BiocCheck RPA_1.17.13.tar.gz 
