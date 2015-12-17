## (Finished,not finished) denoted by (+,-)
# -01. Give more valid starting points
# -02. Check figures

#############################################################################
################## List of log-concave distributions ########################
## 0. Built-in distributions in R (dnorm, dgamma, dbeta, dexp, )
## 1. Normal 
## 2. Gamma if shape parameter>=1
## 3. Beta if both shape parameters>=1
## 4. Exponential 
## 5. Chi Square if df>=2
## 6. Logistic 
## 7. Uniform
## 8. Weibull if shape parameter>=1
## 9. Double Exponential (Laplace)
## 10. Wishart, Dichlet, extreme value distribution

## The following distributions are non-log-concave for all parameters:
## :Pareto, Log normal, Student's t, F-distribution, Cauchy


#############################################################################
########################## 1. Normal distribution ###########################
## 1-1. mean=0 sd=1

### Summary
mysample = ars(100,dnorm,x.start=c(-4,0,4),plot.type="bounds")
mysample = ars(100,dnorm,x.start=c(-4,0,4),plot.type="acceptance")
mysample = ars(100,dnorm,plot.type="bounds")
mysample = ars(100,dnorm,plot.type="acceptance")

### Test: Valid
test_that("sampling ok",{
  
  # Normal distribution
  n <- 100
  x <- rnorm(n,mean=0,sd=1)
  y <- ars(n,dnorm,mean=0,sd=1,x.start=c(-4,0,4))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

### Test: Valid
test_that("sampling ok",{
  
  # Normal distribution
  n <- 100
  x <- rnorm(n,mean=0,sd=1)
  y <- ars(n,dnorm,mean=0,sd=1)
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

########################## 2. Gamma distribution ############################
## 2-1. shape=2 scale=0.5

### Summary
mysample = ars(n = 100,dgamma, shape=2, scale=0.5, x.start=c(0.1,2,5),
               plot.type="bounds")
mysample = ars(n = 100,dgamma, shape=2, scale=0.5, x.start=c(0.1,2,5),
               plot.type="acceptance")
mysample = ars(n = 100,dgamma, shape=2, scale=0.5, plot.type="bounds", domain=c(0,Inf))
mysample = ars(n = 100,dgamma, shape=2, scale=0.5, plot.type="acceptance", domain=c(0,Inf))

### Test without given starting values: Invalid
test_that("sampling ok",{
  x <- rgamma(100,shape=2,scale=0.5)
  y <- ars(100,dgamma,shape=2,scale=0.5,domain=c(0,Inf))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

### Test with given starting values: Invalid
test_that("sampling ok",{
  x <- rgamma(100,shape=2,scale=0.5)
  y <- ars(100,dgamma,shape=2,scale=0.5,x.start=c(0.1,2,5))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

########################## 3. Beta distribution #############################
## 3-1. shape1=1.3 shape2=2.7

### Summary
mysample = ars(n = 100,dbeta, shape1=1.3, shape2=2.7, x.start=c(0.1,0.5,0.9),
               plot.type="bounds")
mysample = ars(n = 100,dbeta, shape1=1.3, shape2=2.7, x.start=c(0.1,0.5,0.9),
               plot.type="acceptance")
mysample = ars(n = 100,dbeta, shape1=1.3, shape2=2.7, plot.type="bounds", domain=c(0,1))
mysample = ars(n = 100,dbeta, shape1=1.3, shape2=2.7, plot.type="acceptance", domain=c(0,1))

### Test without given starting values: Invalid
test_that("sampling ok",{
  x <- rbeta(100,shape1=1.3,shape2=2.7)
  y <- ars(100,dbeta,shape1=1.3,shape2=2.7,domain=c(0,1))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

### Test without given starting values: Invalid
test_that("sampling ok",{
  x <- rbeta(100,shape1=1.3,shape2=2.7)
  y <- ars(100,dbeta,shape1=1.3,shape2=2.7, x.start=c(0.05,0.2,0.95))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

########################## 4. Exponential distribution ######################
## 4.1 rate=2
### Summary
mysample = ars(n = 100,dexp, rate=2, x.start=c(0.1,2,5),plot.type="bounds")
mysample = ars(100,dexp, rate=2, x.start=c(0.1,2,5),plot.type="acceptance")
mysample = ars(n = 100, dexp, rate=2 ,domain=c(0,Inf), plot.type="bounds")
mysample = ars(100, dexp, rate=2 ,domain=c(0,Inf), plot.type="acceptance")

### Test without given starting values: Valid
test_that("sampling ok",{
  x <- rexp(100,rate=2)
  y <- ars(100,dexp,rate=2,domain=c(0,Inf))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

### Test with given starting values: Valid
test_that("sampling ok",{
  x <- rexp(100,rate=2)
  y <- ars(100,dexp, rate=2,x.start=c(0.1,1,2))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})


############################## 5. Chi #######################################
## 5.1 df=2
### Summary - Right plot?
mysample = ars(100,dchisq,df=2,x.start=c(0.1,2,4),plot.type="bounds")
mysample = ars(100,dchisq,df=2,x.start=c(0.1,2,4),plot.type="acceptance")
mysample = ars(100,dchisq,df=2,domain=c(0,Inf), plot.type="bounds")
mysample = ars(100,dchisq,df=2,domain=c(0,Inf), plot.type="acceptance")

### Test without given starting values: Invalid
test_that("sampling ok",{
  x <- rchisq(100,df=2)
  y <- ars(100,dchisq,df=2,domain=c(0,Inf))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

### Test with given starting values: Valid
test_that("sampling ok",{
  x <- rchisq(100,df=2)
  y <- ars(100,dchisq,df=2,x.start=c(0.1,2,4))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})


############################## 6. Logistic ###################################
## 6.1 location=0 scale=1
### Summary
mysample = ars(n = 100,dlogis, location=0, scale=1, x.start=c(-2,0,2),plot.type="bounds")
mysample = ars(100,dlogis, location=0, scale=1, x.start=c(-2,0,2),plot.type="acceptance")
mysample = ars(n = 100, dlogis, location=0, scale=1, plot.type="bounds")
mysample = ars(100, dlogis, location=0, scale=1, plot.type="acceptance")

### Test without given starting values: Valid
test_that("sampling ok",{
  x <- rlogis(100, location=0, scale=1)
  y <- ars(100,dlogis,location=0, scale=1)
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

### Test with given starting values: Valid
test_that("sampling ok",{
  x <- rlogis(100, location=0, scale=1)
  y <- ars(100,dlogis,location=0, scale=1,x.start=c(-2,0,2))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

############################## 7. Uniform ###################################
## 7.1 min=0 max=1
### Summary
mysample = ars(n = 100,dunif, min=0, max=1, x.start=c(0.1,0.5,0.9),plot.type="bounds")
mysample = ars(100,dunif, min=0, max=1, x.start=c(0.1,0.5,0.9),plot.type="acceptance")
mysample = ars(n = 100, dunif, domain=c(0,1), plot.type="bounds")
mysample = ars(100, dunif, domain=c(0,1), plot.type="acceptance")

### Test without given starting values: Valid
test_that("sampling ok",{
  x <- runif(100)
  y <- ars(100,dunif,domain=c(0,1))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

### Test with given starting values: Valid
test_that("sampling ok",{
  x <- runif(100)
  y <- ars(100,dunif,x.start=c(0.1,0.5,0.9))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

################## 8. Weibull if shape parameter>=1 ###################################
## 8.1 shape=1 scale=1
### Summary
mysample = ars(n = 100,dweibull, shape=1, scale=1, x.start=c(0.001,2,4),plot.type="bounds")
mysample = ars(100,dweibull, shape=1, scale=1, x.start=c(0.001,2,4),plot.type="acceptance")
mysample = ars(n = 100, dweibull, shape=1, domain=c(0,Inf), plot.type="bounds")
mysample = ars(100, dweibull, shape=1, domain=c(0,Inf), plot.type="acceptance")

### Test without given starting values: Invalid
test_that("sampling ok",{
  x <- rweibull(100,shape=1)
  y <- ars(100,dweibull,shape=1,domain=c(0,Inf))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

### Test with given starting values: Valid
test_that("sampling ok",{
  x <- rweibull(100,shape=1)
  y <- ars(100,dweibull,shape=1,x.start=c(0.001,2,4))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})



############## 9. Double Exponential (Laplace) #######################################
## 9.1 mu=0 lambda=1
install.packages("smoothmest")
library(smoothmest)

### Summary
mysample = ars(n = 100,ddoublex, mu=0, lambda=1, x.start=c(-2,0,2),plot.type="bounds")
mysample = ars(100,ddoublex, mu=0, lambda=1, x.start=c(-2,0,2),plot.type="acceptance")
mysample = ars(n = 100,ddoublex, mu=0, lambda=1, domain=c(-Inf,Inf), plot.type="bounds")
mysample = ars(100,ddoublex, mu=0, lambda=1, domain=c(-Inf,Inf), plot.type="acceptance")

### Test without given starting values: Valid
test_that("sampling ok",{
  x <- rdoublex(100)
  y <- ars(100,ddoublex,domain=c(-Inf,Inf))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

### Test with given starting values: Valid
test_that("sampling ok",{
  x <- rdoublex(100)
  y <- ars(100,ddoublex,x.start=c(-2,0,2))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})



########################### Other examples ############################################
# Example 01: concave parabola, valid
g = function(x) {exp(-(1/2)*(x)^2)}
test.sample = ars(100,g)
test.sample = ars(100,g,x.start=c(-2,0,2))
test.sample = ars(100,g,x.start=c(-2,0,2),plot.type="bounds")
test.sample = ars(100,g,x.start=c(-2,0,2),plot.type="acceptance")
test.sample = ars(100,g,plot.type="acceptance")

# Example 04: sin(x) in 1:6, non valid
# mysample = ars(n=100,g = sin, dg = cos,x.start=c(3,1,5),plot.type="bounds")
# mysample = ars(n=100,g = sin, x.start=c(3,1,5),plot.type="acceptance")

# Example 05: sin(x) in 1:3, valid
mysample = ars(n=100, g = sin, dg = cos,x.start=c(3,1,2,0.1),plot.type="bounds")
mysample = ars(n=100,g = sin, x.start=c(3,1,2,0.1),plot.type="acceptance")
## mysample = ars(n=100, g = sin, dg = cos, plot.type="bounds") 
## mysample = ars(n=100,g = sin, plot.type="acceptance")


########################### Efficiency check #######################################
## Comparing computation speeds between our function and the built-in ars
g = function(x) {-0.5*x^2}
dg = function(x) {-x}
## Our function
system.time(ars(100, g=g,dg=dg,log.transform=TRUE,x.start=c(-4,1,4)))
# user  system elapsed 
# 0.00    0.00    0.02
rm(list=ls())
g = function(x) {-0.5*x^2}
dg = function(x) {-x}
## Built-in ARS
system.time(ars(1000,f=g,fprima=dg))
# user  system elapsed 
# 0.000   0.000   0.163 
