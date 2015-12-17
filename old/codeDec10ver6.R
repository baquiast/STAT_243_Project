###################
## Going further ##
###################

# 1. Option to ars() could be log.transformed.f = F
# 2. Derivative of the function can be given by users as: => DONE
# 3. Figure out what to do when x.cand is element in x.vec => DONE
# 4. Give the option to give function as input in another way 
# 5. Find a way to do the convexity check of f by evaluating fewer f
# 6. Take unbounded domains into consideration => DONE in ver7: ars() in R has starting points as input, can we do in the same way?
# 7. Include the convexity check => DONE
# 8. Can we set the three starting points as a user input and default of the points as (-4,1,4) (to make our life easier)? => DONE in ver7

##############
## Comments ##
##############

# 1. "f" is the log of the users density function
# 2. "dom" is the domain defined for the density
# 3. "n" is the requested number of sample elements
# 4. The candidates, x.cand, for a sample value are drawn from the density
#     defined by the upper bound of f. First x.cand is drawn from a equidistant
#     set in dom, with steps 2*dx. The candidate is then adjusted by U(-dx,dx),
#     in order to cover dom.
# 5. Setting plot.type = "acceptance" depend on evaluating f a lot of times and
#     should not be done when f is expensive to evaluate.

###################
## Sub-functions ##
###################
library(numDeriv)

# f convexity check
is.convex = function(f,dom){
  x.seq <- seq(dom[1],dom[2], l=100)
  for (i in 3:100){
    if((f(x.seq[i])-f(x.seq[i-1])) > (f(x.seq[i-1])-f(x.seq[i-2]))) stop("Error: The distribution is not log concave!")
  }
}

# Evaluate lower bound at x
lower = function(x, x.vec, f.vec){
  l = max(which(x>x.vec))
  return(f(x.vec[l])+(x-x.vec[l])*((f.vec[l+1]-f.vec[l])/(x.vec[l+1]-x.vec[l])))
}

# Evaluate upper bound at x
upper = function(x, x.vec, f.vec, d.vec){
  l = max(which(x>x.vec))
  return(min(f.vec[l]+(x-x.vec[l])*d.vec[l],f.vec[l+1]+(x-x.vec[l+1])*d.vec[l+1]))
}

# Plot of the acceptance region and bounds
boundary.plot = function(x.vec,f.vec,d.vec,x.temp,with.f = F, f=F){
  # Create upper bound
  nvec = length(x.vec)
  x.center = rep(NA,nvec-1)
  f.center = x.center
  for (i in 1:(nvec-1)){
    a=sapply(x.temp, function(x){
      f.vec[i]+(x-x.vec[i])*d.vec[i]>(f.vec[i+1]+(x-x.vec[i+1])*d.vec[i+1])})
    x.center[i] = x.temp[which(diff(a)!=0)]
    f.center[i] = f.vec[i]+(x.center[i]-x.vec[i])*d.vec[i]
  }
  x.corners = c(x.vec[1],x.center,x.vec[nvec])
  f.corners = c(f.vec[1],f.center,f.vec[nvec])
  
  # Plot
  if(with.f==T){
    plot(f, xlim=dom, lwd =1)
    polygon(c(x.corners,x.vec[1]),c(f.corners,f.vec[1]),col="firebrick")
    polygon(c(x.temp,x.vec[1]),c(f(x.temp),f.vec[1]),col="palegreen4")
    lines(x.corners,f.corners, lwd=1)
  } else{
    plot(x.corners,f.corners, type = "l",lty =1, lwd=1, col="red")
  }
  lines(x.vec,f.vec, lwd =1)
}

###################
## Main function ##
###################

ars = function(f, df = NA, x.start = NA, dom, n, plot.type = "none"){ 
  
  set.seed(200)
  is.convex(f,dom) 
  x.out=rep(NA,n)
  if(length(x.start) == 3) {x.vec = x.start} else {x.vec = c(dom[1],mean(c(dom[1],dom[2])),dom[2])}
  f.vec = f(x.vec)
  if(class(df)=="function") {d.vec=df(x.vec)} else {d.vec = grad(f,x.vec)}
  dx = 1e-2
  x.temp = seq(x.vec[1]+dx,x.vec[3]-dx,2*dx)
  i = 1
  
  while(i<=n){
    u = runif(1,0,1)
    x.cand = sample(x.temp[!x.temp %in% x.out],size=1,prob=sapply(x.temp[!x.temp %in% x.out], function(x){
      exp(upper(x,x.vec,f.vec,d.vec))}))
    x.cand = runif(1,x.cand-dx,x.cand+dx)
    
    if(u <= exp(lower(x.cand,x.vec,f.vec)-upper(x.cand,x.vec,f.vec,d.vec))) {
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
      if(class(df)=="function") {d.vec = append(d.vec,df(x.cand),after=l)} else {d.vec = append(d.vec,grad(f,x.cand),after=l)}
    }
  }
  
  if(plot.type == "bounds"){
    boundary.plot(x.vec,f.vec,d.vec, x.temp, with.f=F)
  }
  if(plot.type == "acceptance"){
    boundary.plot(x.vec,f.vec,d.vec, x.temp, with.f=T, f)
  }
  
  return(x.out)
}


##########
## test ##
##########

# User input
f = function(x) {-(1/2)*(x)^2}
df = function(x) {-x}
dom = c(-10,10)
n = 100

# Draw a sample using ars()
test.sample = ars(f,df=NA,x.start=NA,dom,n)
test.sample = ars(f,df=NA,x.start=NA,dom,n,plot.type="bounds")
test.sample = ars(f,df=NA,x.start=NA,dom,n,plot.type="acceptance")

test.sample = ars(f,df,x.start=c(-4,1,4),dom,n) # Domain and starting points are repetitive information?

hist(test.sample)


#########################################
## Tests for complicated distributions ##
#########################################

# Example 01: gamma(2,0.5)
f <- function(x,shape=2,scale=0.5){(shape-1)*log(x)-x/scale}
df <- function(x,shape=2,scale=0.5) {(shape-1)/x-1/scale}
dom = c(0.001,10)
n = 100
mysample <- ars(f,df,x.start=NA,dom,n)
mysample
hist(mysample) 
hist(rgamma(100,shape=2,scale=0.5))

#Example 02: beta(1.3,2.7)
f <- function(x,a=1.3,b=2.7){(a-1)*log(x)+(b-1)*log(1-x)}
df <- function(x,a=1.3,b=2.7){(a-1)/x-(b-1)/(1-x)}
dom = c(0.001,0.999)
mysample <- ars(f,df,x.start=NA,dom,n)
mysample
hist(mysample)
hist(rbeta(100,shape1=1.3,shape2=2.7))
