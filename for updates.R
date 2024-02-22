setwd(getwd())
setwd("C:/Users/Lenovo/Documents/GitHub")
library(Rcpp)

#run to build the skeleton of the package - but maybe not serve anymore
#Rcpp.package.skeleton("")

#run when the cpp (and maybe also R?) function is modified
Rcpp::compileAttributes("C:/Users/Lenovo/Documents/GitHub/GiniDecA")

#adjust DESCRIPTION and NAMESPACE to import necessary packages

devtools::install_github("FedericoAttili/GiniDecA",force = T)

library(GiniDecA)
giniDec(c(1,1,2,2),c(1,1,2,2))
