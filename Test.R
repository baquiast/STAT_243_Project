# Testing sample ok

library(testthat)

n <- 1000

# Code for testing the sample
test_sample <- function(x,y){
  test_passed <- TRUE
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
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



# Normal distribution
test_that("Normal Distribution",{
  x <- rnorm(n,mean=1,sd=1)
  y <- ars(n,dnorm,mean=1,sd=1,x.start=c(-10,0,10))
  test_passed <- test_sample(x,y)
  expect_true(test_passed)
  if(test_passed) {print('Normal distribution well sampled')}
})

# Beta distribution with alpha = beta = 2
test_that("Beta distribution",{
  x <- rnorm(n,mean=1,sd=1)
  y <- ars(n,dbeta,shape1=2,shape2=2,domain=c(0,1))
  test_passed <- test_sample(x,y)
  expect_true(test_passed)
  if(test_passed) {print('Beta distribution well sampled')}
})
