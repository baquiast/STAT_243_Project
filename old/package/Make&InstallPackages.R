## Make pacakges
devtools::load_all()
devtools::document()
R CMD build ars

## Install packages
install.packages("/home/oski/R/ars_0.1.2.tar.gz", repos = NULL, type="source")
library(ars)

mysample = ars(100, dnorm, x.start=c(-4,0,4), plot.type="bounds")
mysample = ars(100, dnorm, x.start=c(-4,0,4), plot.type="acceptance")
mysample = ars(100,dnorm,plot.type="bounds")
mysample = ars(100,dnorm,plot.type="acceptance")

library(testthat)
test_package('ars','main')
library(RUnit)
test_package('ars','units')
