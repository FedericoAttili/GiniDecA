setwd(getwd())
library(Rcpp)

#run to build the skeleton of the package - but maybe not serve anymore
#Rcpp.package.skeleton("")

#run when the cpp (and maybe also R?) function is modified
Rcpp::compileAttributes(getwd())


#adjust DESCRIPTION and NAMESPACE to import necessary packages
devtools::document()

devtools::install_github("FedericoAttili/GiniDecA",force = T)

library(GiniDecA)
giniDec(c(1,1,2,2),c(1,1,2,2))
