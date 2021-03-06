###################
## Going further ##
###################

## (Finished,not finished) denoted by (+,-)
# +01. Option to ars() could be log.transformed.f = F
# +02. Derivative of the function can be given by users as
# +03. Figure out what to do when x.cand is element in x.vec
# +04. Give the option to give function as input in another way 
# +05. Find a way to do the convexity check of f by evaluating fewer f
# -06. Take unbounded domains into consideration => in ver7: ars() in 
#       R has starting points as input, can we do in the same way?
# -07. Give error when starting point is not compatible with density
# +08. Give error when we cannot take log of output from density in specified
#     interval
# +09. Give error when starting values are not (unique, numbers,>=3)
# -10. Figure out if/how the user should have the necessary packages installed 

##############
## Comments ##
##############

# 01. "f" is the log of the users density function
# 02. "dom" is the domain defined for the density
# 03. "n" is the requested number of sample elements
# 04. The candidates, x.cand, for a sample value are drawn from the density
#     defined by the upper bound of f. First x.cand is drawn from a equidistant
#     set in dom, with steps 2*dx. The candidate is then adjusted by U(-dx,dx),
#     in order to cover dom.
# 05. Setting plot.type = "acceptance" depend on evaluating f a lot of times and
#     should not be done when f is expensive to evaluate.
# 06. x.start is a vector of starting points
# 07. 

###################
## Sub-functions ##
###################
library(numDeriv)
library(testthat)

# Check for concavity across three (x,log(f(x)))
is.local.concave = function(f.vals,x.vals){
  if(((f.vals[3]-f.vals[2])/(x.vals[3]-x.vals[2]))>
       ((f.vals[2]-f.vals[1])/(x.vals[2]-x.vals[1]))){
    stop("Error: The distribution is not log concave!")
  }
}

# Evaluate lower bound at x
lower = function(x, x.vec, f.vec){
  l = max(which(x>x.vec))
  return(f.vec[l]+(x-x.vec[l])*((f.vec[l+1]-f.vec[l])/(x.vec[l+1]-x.vec[l])))
}

# Evaluate upper bound at x
upper = function(x, x.vec, f.vec, d.vec){
  l = max(which(x>x.vec))
  return(min(f.vec[l]+(x-x.vec[l])*d.vec[l],
             f.vec[l+1]+(x-x.vec[l+1])*d.vec[l+1]))
}

# Plot of the acceptance region and bounds
boundary.plot = function(x.vec,f.vec,d.vec, f=FALSE){
  # Create upper bound
  nvec = length(x.vec)
  x.vals = seq(x.vec[1], x.vec[nvec], l = 1000)
  x.center = rep(NA,nvec-1)
  f.center = x.center
  for (i in 1:(nvec-1)){
    a=sapply(x.vals, function(x){
      f.vec[i]+(x-x.vec[i])*d.vec[i]>(f.vec[i+1]+(x-x.vec[i+1])*d.vec[i+1])})
    idx = ifelse(sum(which(diff(a)!=0))==0,length(x.vals)/2,which(diff(a)!=0))
    x.center[i] = x.vals[idx]
    f.center[i] = f.vec[i]+(x.center[i]-x.vec[i])*d.vec[i]
  }
  x.corners = c(x.vec[1],x.center,x.vec[nvec])
  f.corners = c(f.vec[1],f.center,f.vec[nvec])
  
  # Plot
  if(is.function(f)){
    plot(f, xlim=c(x.vec[1],x.vec[nvec]), lwd =1, yaxt='n', ylab='', xlab="x",
         main ="Acceptance region and rejection of the log density",
         cex.lab = 1.4)
    polygon(c(x.corners,x.vec[1]),c(f.corners,min(f.vec[1],f.vec[nvec])),
            col="firebrick", border="firebrick")
    polygon(c(x.vals,x.vec[nvec],x.vec[1]),c(f(x.vals),min(f.vec),min(f.vec)),
            col="palegreen4",border="palegreen4")
    lines(x.corners,f.corners, lwd=1, col="firebrick")
    lines(x.vec,f.vec, lwd =1, col="grey10")
    legend(x="topright", legend=c("Acceptance region", "Rejection region",
                      "Lower bound"),col=c("palegreen4","firebrick","grey10"),
           pch=c(15,15,NA),lty=c(0,0,1), seg.len = 0.5, cex = 0.6)
  } else{
    plot(x.corners,f.corners, type = "l",lty =1, lwd=1, col="red", xlab="x",
         yaxt='n', ylab='', main = "Upper and lower bound of the log density",
         cex.lab = 1.4)
    lines(x.vec,f.vec, lwd =1,col="black")
    legend(x="topright",legend=c("Upper bound", "Lower bound"),
           col=c("black","red"), lty=c(1,1), cex = 0.6)
  }
}

###################
## Main function ##
###################

ars = function(n, g, ..., dg=NULL, x.start=c(-4,1,4), plot.type = "none"){ 
  
  # Check if there is at least three starting values
  if(length(x.start)<3){
    stop("At least three starting values are required")
  }
  
  # Check if starting values are unique
  if(length(unique(x.start))!=length(x.start)){
    stop("The starting values are not unique")
  }
  
  # Log transformations
  f = function(x){
    if(sum(is.na(log(g(x,...))))!=0){
      stop("The input density has negative values")
    }
    return(log(g(x,...)))
  }
  
  df = function(x){
    dg(x)/g(x)
  }
  
  # Initialization
  set.seed(200)
  x.out=rep(NA,n)
  x.vec = sort(x.start)
  f.vec = f(x.vec)
  if(is.function(dg)){d.vec=df(x.vec)} else {d.vec = grad(f,x.vec)}
  dx = 1e-2
  x.temp = seq(x.vec[1]+dx,x.vec[length(x.vec)]-dx,2*dx)
  i = 1
  
  # Create vector of samples
  while(i<=n){
    u = runif(1,0,1)
    x.cand = sample(x.temp,size=1,prob=sapply(x.temp, function(x){
      exp(upper(x,x.vec,f.vec,d.vec))}))
    x.cand = runif(1,x.cand-dx,x.cand+dx)
    
    if(u <= exp(lower(x.cand,x.vec,f.vec)-upper(x.cand,x.vec,f.vec,d.vec))){
      x.out[i] = x.cand
      i=i+1
    }else{
      f.cand = f(x.cand)
      if(u <= exp(f.cand-upper(x.cand,x.vec,f.vec,d.vec))){
        x.out[i] = x.cand
        i = i+1
      }
      l = max(which(x.cand>x.vec))
      x.vec = append(x.vec,x.cand,after=l)
      f.vec = append(f.vec,f.cand,after=l)
      if(is.function(dg)){d.vec = append(d.vec,df(x.cand),after=l)}else{
        d.vec = append(d.vec,grad(f,x.cand),after=l)}
      is.local.concave(f.vec[l:(l+2)],x.vec[l:(l+2)])
    }
  }
  
  # Optional plot part
  if(plot.type == "bounds"){
    boundary.plot(x.vec,f.vec,d.vec)
  }
  if(plot.type == "acceptance"){
    boundary.plot(x.vec,f.vec,d.vec,f)
  }
  
  return(x.out)
}

###########################
## Visual Representation ##
###########################

n <- 1000
x <- rnorm(n,mean=0,sd=1)
y <- ars(n, dnorm, mean=0,sd=1, df=NA, x.start=c(-10,1,10), plot.type="bounds")
x0 <- c(min(x,y),max(x,y))
qqplot(x, y, plot.it = TRUE)
lines(x0,x0,col="blue")

###############
## Test that ##
###############

test_that("sampling ok",{
  
  # Normal distribution
  n <- 1000
  x <- rnorm(n,mean=1,sd=1)
  y <- ars(n,dnorm,mean=1,sd=1,x.start=c(-10,0,10))
  wilcox_test <- wilcox.test(x, y, 
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
                             conf.int = FALSE, conf.level = 0.95)$p.value
  ks <- ks.test(x, y, alternative = "two.sided",exact = NULL)$p.value
  expect_true(wilcox_test>0.05)
  expect_true(ks>0.05)
})

# Test loop within test_that (Guillaume)
test <- c(parse(text='1+1'))
eval(test[1])

#distributions <- c(parse(text='log_dnorm'))
#ars(10,eval(distributions[1]),x.start=c(-10,0,10))

##########
## test ##
##########

# Example 01: concave parabola, valid
n = 100
f = function(x) {exp(-(1/2)*(x)^2)}
test.sample = ars(n,f)
test.sample = ars(n,f,x.start=c(-2,0,2))
test.sample = ars(n,f,x.start=c(-2,0,2),plot.type="bounds")
test.sample = ars(n,f,x.start=c(-2,0,2),plot.type="acceptance")

# Example 02: gamma(2,0.5), valid
mysample = ars(n = 100,dgamma, shape=2, scale=0.5, x.start=c(0.1,2,5),
               plot.type="bounds")
mysample = ars(n,dgamma, shape=2, scale=0.5, x.start=c(0.1,2,5),
               plot.type="acceptance")

# Example 03: beta(1.3,2.7), valid
mysample = ars(n = 100, dbeta, shape1 = 1.3, shape2 = 2.7, 
               x.start=c(0.1,0.5,0.9), plot.type="bounds")
mysample = ars(n, dbeta, shape1 = 1.3, shape2 = 2.7, 
               x.start=c(0.1,0.5,0.9), plot.type="acceptance")

# Example 04: sin(x) in 1:6, non valid
# mysample = ars(n=100,g = sin, dg = cos,x.start=c(3,1,5),plot.type="bounds")
# mysample = ars(n=100,g = sin, x.start=c(3,1,5),plot.type="acceptance")

# Example 05: sin(x) in 1:3, valid
mysample = ars(n=100, g = sin, dg = cos,x.start=c(3,1,2,0.1),plot.type="bounds")
mysample = ars(n=100,g = sin, x.start=c(3,1,2,0.1),plot.type="acceptance")

# Example 06: dnorm(x) in -2,2, valid
mysample = ars(100,dnorm,x.start=c(1,-2,3,5),plot.type="bounds")
mysample = ars(100,dnorm,x.start=c(1,-2,3,5),plot.type="acceptance")

# Example 07: exp(2), valid
mysample = ars(n = 100,dexp, rate=2, x.start=c(0.1,2,5),plot.type="bounds")
mysample = ars(n,dexp, rate=2, x.start=c(0.1,2,5),plot.type="acceptance")

