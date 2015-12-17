library(numDeriv)

test.is.local.concave <- function() {
  h_Normal = function(x) {log(dnorm(x,0,1))}
  x_Normal <- c(-4,0,4)
  d_Normal <- grad(h_Normal,x_Normal)
  checkEquals(is.local.concave(d_Normal),"TRUE")
  
  h_Exp = function(x) {log(dexp(x,1))}
  x_Exp <- c(1,4,100)
  d_Exp <- grad(h_Exp,x_Exp)
  checkEquals(is.local.concave(d_Exp),"TRUE")
  
  h_ChiSq = function(x) {log(dchisq(x,2))}
  x_ChiSq <- c(1,4,10)
  d_ChiSq <- grad(h_ChiSq,x_ChiSq)
  checkEquals(is.local.concave(d_ChiSq),"TRUE")
  
  h_Beta = function(x) {log(dbeta(x,3,3))}
  x_Beta <- c(0.2,0.4,0.9)
  d_Beta <- grad(h_Beta,x_Beta)
  checkEquals(is.local.concave(d_Beta),"TRUE")
  
  h_t = function(x) {log(dt(x,2))}
  x_t <- c(0,1,5)
  d_t <- grad(h_t,x_t)
  checkException(is.local.concave(d_t), "Error: The distribution is not log concave!")
  
  h_Cauchy = function(x) {log(dcauchy(x,1,0.5))}
  x_Cauchy <- c(0,0.2,2)
  d_Cauchy <- grad(h_Cauchy,x_Cauchy)
  checkException(is.local.concave(d_Cauchy), "Error: The distribution is not log concave!")
}

test.is.bounded <- function() {
  x.vec <- c(-4,-1.091335,0,1.548574,4)
  f.vec <- c(8.9189385,-1.5144445,-0.9189385,-2.1179792,-8.9189385)
  d.vec <- c(4,1.091335,0,-1.548574,-4)
  x.cand.re0 <- 2
  f.cand.re0 <- log(dnorm(2))
  checkEquals(is.bounded(x.cand.re0,x.vec,f.vec,d.vec,f.cand.re0),"TRUE")
}