test.gen.start <- function() {
  x.start=0
  domain=c(-4,4)
  f=function(x) {dnorm(x,0,1)}
  checkException(gen.start(x.start,domain,f), "At least three starting values are required")
  x.start=c(4,0,4)
  checkException(gen.start(x.start,domain,f), "The starting values are not unique")
  x.start=c()
  domain=0
  checkException(gen.start(x.start,domain,f), "Invalid domain")
  domain=c(4,4)
  checkException(gen.start(x.start,domain,f), "Invalid domain")
  domain=c(-4,0,4)
  checkException(gen.start(x.start,domain,f), "Invalid domain")
  domain=c(-1,1)
  checkEquals(gen.start(x.start,domain,f),c(-1,-8.326673e-17,1),tolerance=10^-6)
  domain=c(-Inf,Inf)
  checkEquals(gen.start(x.start,domain,f), c(-1e+05,1e+05,1e+05),tolerance=10^-6)
  f=function(x) {dexp(x,1)}
  domain=c(0,10)
  checkEquals(gen.start(x.start,domain,f), c(0,5.575865e-05,10),tolerance=10^-6)
  domain=c(0,Inf)
  checkEquals(gen.start(x.start,domain,f), c(0e+00,1e+05,1e+05),tolerance=10^-6)
  f=function(x) {dchisq(x,2)}
  domain=c(0,10)
  checkEquals(gen.start(x.start,domain,f), c(0,5.575865e-05,10),tolerance=10^-6)
  domain=c(0,Inf)
  checkEquals(gen.start(x.start,domain,f), c(0e+00,1e+05,1e+05),tolerance=10^-6)
}