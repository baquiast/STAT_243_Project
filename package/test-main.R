############################################################################
## The codes for the main test #############################################
test_sample <- function(x,y){
  test_passed <- TRUE
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, 
                             correct = TRUE, conf.int = FALSE, 
                             conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  
  # Wilcox test
  if(!(wilcox_test>0.05)) {
    test_passed <- FALSE
  }
  
  # KS test
  if(!(ks>0.05)) {
    test_passed <- FALSE
  }
  
  # Final Message
  return(test_passed)
}


#############################################################################
######### List of log-concave distributions to do main tests ################
## 1. Normal 
## 2. Gamma if shape parameter>=1
## 3. Beta if both shape parameters>=1
## 4. Exponential 
## 5. Chi Square if df>=2
## 6. Logistic 
## 7. Uniform
## 8. Weibull if shape parameter>=1
## 9. Double Exponential (Laplace)



## 10. Extreme value distribution

## The following distributions are non-log-concave for all parameters:
## :Pareto, Log normal, Student's t, F-distribution, Cauchy


#############################################################################
## 1. Normal distribution(Mean = 0 SD = 1) ##################################

test_that("Normal Distribution",{
  set.seed(100)
  x = rnorm(n=100,mean=0,sd=1)
  y = ars(n=100,dnorm,mean=0,sd=1,x.start=c(-4,4), plot.type="acceptance")
  test_passed = test_sample(x,y)
  expect_true(test_passed)
  if(test_passed) {
    print('Normal distribution: The main test is passed.')
  }
})

#############################################################################
## 2. Gamma distribution(shape = 1, scale = 1) ##############################

test_that("Gamma Distribution",{
  set.seed(100)
  x = rgamma(n=100,shape=1,scale=1)
  y = ars(n=100,dgamma,shape=1,scale=1,domain=c(0,Inf))
  test_passed = test_sample(x,y)
  expect_true(test_passed)
  if(test_passed) {
    print('Gamma distribution: The main test is passed.')
  }
})

#############################################################################
## 3. Beta distribution(shape1 = 1, shape2 = 1) #############################

test_that("Beta Distribution",{
  set.seed(100)
  x = rbeta(n=100,shape1=1,shape2=1)
  y = ars(n=100,dbeta,shape1=1,shape2=1,domain=c(0,1),plot.type ="acceptance")
  test_passed = test_sample(x,y)
  expect_true(test_passed)
  if(test_passed) {
    print('Beta distribution: The main test is passed.')
  }
})

#############################################################################
## 4. Exponential distribution(rate = 2) ####################################

test_that("Exponential Distribution",{
  set.seed(100)
  x = rexp(n=100,rate=2)
  y = ars(n=100,dexp,rate=2,domain=c(0,Inf))
  test_passed = test_sample(x,y)
  expect_true(test_passed)
  if(test_passed) {
    print('Exponential distribution: The main test is passed.')
  }
})

#############################################################################
## 5. Chi Square distribution(df = 2) #######################################

test_that("Chi Square Distribution",{
  set.seed(100)
  x = rchisq(n=100,df=2)
  y = ars(n=100,dchisq,df=2,domain=c(0,Inf))
  test_passed = test_sample(x,y)
  expect_true(test_passed)
  if(test_passed) {
    print('Chi square distribution: The main test is passed.')
  }
})

#############################################################################
## 6. Logistic distribution #################################################

test_that("Logistic Distribution",{
  set.seed(100)
  x = rlogis(n=100, location=0, scale=1)
  y = ars(n=100, f=dlogis, location=0, scale=1)
  test_passed = test_sample(x,y)
  expect_true(test_passed)
  if(test_passed) {
    print('Logistic distribution: The main test is passed.')
  }
})

#############################################################################
## 7. Uniform distribution ##################################################

test_that("Uniform Distribution",{
  set.seed(100)
  x = runif(n=100, min=0, max=1)
  y = ars(n=100, f=dunif, min=0, max=1, domain=c(0,1))
  test_passed = test_sample(x,y)
  expect_true(test_passed)
  if(test_passed) {
    print('Uniform distribution: The main test is passed.')
  }
})

#############################################################################
## 8. Weibull distribution ##################################################

test_that("Weibull Distribution",{
  set.seed(100)
  x = rweibull(n=100, shape=1, scale=1)
  y = ars(n=100, f=dweibull, shape=1, scale=1, domain=c(0,Inf))
  test_passed = test_sample(x,y)
  expect_true(test_passed)
  if(test_passed) {
    print('Weibull distribution: The main test is passed.')
  }
})

#############################################################################
## 9. Double Exponential distribution #######################################

library(smoothmest)
test_that("Double Exponential Distribution",{
  set.seed(100)
  x = rdoublex(n=100,  mu=0, lambda=1)
  y = ars(n=100, f=ddoublex, mu=0, lambda=1)
  test_passed = test_sample(x,y)
  expect_true(test_passed)
  if(test_passed) {
    print('Double Exponenetial distribution: The main test is passed.')
  }
})
