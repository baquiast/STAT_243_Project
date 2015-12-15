test.gen.start <- function() {
  x.start=0
  domain=c(-4,4)
  f=function(x) {log(dnorm(x,0,1))}
  checkException(gen.start(x.start,domain,f), "At least three starting values are required")
  x.start=c()
  domain=0
  checkException(gen.start(x.start,domain,f), "Invalid domain")
  domain=c(4,4)
  checkException(gen.start(x.start,domain,f), "Invalid domain")
  domain=c(-4,0,4)
  checkException(gen.start(x.start,domain,f), "Invalid domain")
  domain=c(-1,1)
  checkEquals(gen.start(x.start,domain,f),c(-1,1),tolerance=10^-6)
  domain=c(-Inf,Inf)
  checkEquals(gen.start(x.start,domain,f),c(-15,15),tolerance=10^-6)
  f=function(x) {log(dexp(x,1))}
  domain=c(0,10)
  checkEquals(gen.start(x.start,domain,f), c(0.000525368,5.000262758),tolerance=10^-6)
  domain=c(0,Inf)
  checkEquals(gen.start(x.start,domain,f), c(0.000525368,7.500262796),tolerance=10^-6)
  f=function(x) {log(dchisq(x,2))}
  domain=c(0,10)
  checkEquals(gen.start(x.start,domain,f), c(0.000525368,5.000262758),tolerance=10^-6)
  domain=c(0,Inf)
  checkEquals(gen.start(x.start,domain,f), c(0.000525368,7.500262796),tolerance=10^-6)
}