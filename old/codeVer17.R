rm(list=ls())

###################
## Going further ##
###################

## (Finished,not finished) denoted by (+,-)
# +01. Option to ars() could be log.transformed.f = F
# +02. Derivative of the function can be given by users
# +03. Figure out what to do when x.cand is element in x.vec
# +04. Give the option to give function as input in another way 
# +05. Find a way to do the convexity check of f by evaluating fewer f
# +06. Take unbounded domains into consideration
# +07. Give error when starting point is not compatible with density
# +08. Give error when we cannot take log of output from density in specified
#     interval
# +09. Give error when starting values are not (unique, numbers,>=3)
# -10. Figure out if/how the user should have the necessary packages installed 
# -11. Is it possible to have one or two starting points? (See built-in ars)
# +12. Vectorize as much as possible
# -13. Efficiency problem

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

# Boundedness
is.bounded = function(x, x.vec, f.vec, d.vec, f){
  if(upper(x, x.vec, f.vec, d.vec) < f(x) | f(x) < lower(x, x.vec, f.vec)){
    stop("Error: The distribution is not bounded!")
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

# Generate starting points
gen.start = function(x.start,domain) {
  # When user input starting values
  if(length(x.start)!=0) {
    # Check if there is at least three starting values
    if(length(x.start)<3){
      stop("At least three starting values are required")
    }
    # Check if starting values are unique
    if(length(unique(x.start))!=length(x.start)){
      stop("The starting values are not unique")
    }
  }else{ # User does not give starting values
    if(length(domain)!=2 | length(unique(domain))!=2){ # Invalid domain
     stop("Invalid domain")
    } else{ # Valid domain
      if(domain[1]>domain[2]){
        domain = c(domain[2],domain[1])
      }
      if(domain[1]==-Inf){
        domain[1] = -1e5
      }
      if(domain[2]==Inf){
        domain[2] = 1e5
      }
    
      # Create starting values from domain
      p1 = seq(0.1,0,length.out=10)
      p2 = seq(10,0,length.out=10)
      p.vec1 = cumsum(exp(p1)/sum(exp(p1)))
      p.vec2 = cumsum(exp(p2)/sum(exp(p2)))
      p.vec = unique(sort(c(p.vec1,p.vec2,0.5*p.vec1,0.5*p.vec2)))
      p.vec = p.vec[p.vec!=0 & p.vec!=1]
    
      x.start = domain   
      for(i in 1:2) {
        x.temp1 = rep(NULL,length(p.vec))
        if(f(x.start[i])==Inf){
          x.temp1 = domain[i] - ((-1)^(i))*min(1e-2,(domain[2]-domain[1]))*p.vec
          x.temp1 = x.temp1[abs(f(x.temp1))!=Inf]
          x.start[i] = x.temp1[!is.na(grad(f,x.temp1))][1]
        }
        else if(f(x.start[i])==-Inf){
          x.temp1 = domain[i] - ((-1)^(i))*(domain[2]-domain[1])*p.vec
          x.temp1 = x.temp1[abs(f(x.temp1))!=Inf]
          x.start[i] = x.temp1[!is.na(grad(f,x.temp1))][1]
        }
        else {
          if(!is.na(grad(f,x.start[i]))) {
            x.start[i] = x.start[i] 
          } else {
            x.temp1 = domain[i] - ((-1)^(i))*(domain[2]-domain[1])*p.vec
            x.start[i] = x.temp1[!is.na(grad(f,x.temp1))][1]            
          }
        }
      }
      x.mod = optimize(f,interval=x.start,maximum=TRUE)$maximum
      x.start = sort(c(x.start, x.mod))
    }
  }
  return(x.start)
}
# Plot of the acceptance region and bounds
boundary.plot = function(x.vec,f.vec,d.vec, f=FALSE){
  # Create upper bound
  nvec = length(x.vec)
  x.vals = seq(x.vec[1], x.vec[nvec], l = 1000)
  idx = f.center = x.center = rep(NA,nvec-1)
  
  mat = cbind(x.vec,f.vec,d.vec)
  idx = sapply(1:(nvec-1), function(i) {
    a = (mat[i,2]+(x.vals-mat[i,1])*mat[i,3])>(mat[i+1,2]+(x.vals-mat[i+1,1])*mat[i+1,3])
    idx = ifelse(sum(which(diff(a)!=0))==0,0,which(diff(a)!=0))
  }
  )
  
  idx[idx==0] =  idx[which(idx==0)+1]
  x.center = x.vals[idx]
  f.center = f.vec[1:(nvec-1)] + (x.center[1:(nvec-1)]-x.vec[1:(nvec-1)]) * d.vec[1:(nvec-1)]
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
            col="palegreen4",border="palegreen4")  # acceptance region
    lines(x.corners,f.corners, lwd=1, col="firebrick")
    lines(x.vec,f.vec, lwd =1, col="grey10") # lower bound
    legend(x="topright", legend=c("Acceptance region", "Rejection region",
                                  "Lower bound"),col=c("palegreen4","firebrick","grey10"),
           pch=c(15,15,NA),lty=c(0,0,1), seg.len = 0.5, cex = 0.6)
  } else{
    plot(x.corners,f.corners, type = "l",lty =1, lwd=1, col="red", xlab="x",
         yaxt='n', ylab='', main = "Upper and lower bound of the log density",
         cex.lab = 1.4)
    lines(x.vec,f.vec, lwd =1,col="black") # lower bound
    legend(x="topright",legend=c("Upper bound", "Lower bound"),
           col=c("black","red"), lty=c(1,1), cex = 0.6)
  }
}

###################
## Main function ##
###################

ars = function(n, g, ..., dg=NULL, log.transform=FALSE, x.start=NULL, domain=c(-Inf,Inf), plot.type = "none"){ 
  
  ## Log transformations 
  if(log.transform==TRUE){
    f = function(x){
      if(sum(is.na(g(x,...)))!=0){
        stop("The input density has negative values or too small values")
      }
      return(g(x,...))
    }  
  } else{
    f = function(x){
      if(sum(is.na(log(g(x,...))))!=0){
        stop("The input density has negative values")
      }
      return(log(g(x,...)))
    }
  }
  
  if(log.transform==TRUE){
    df = function(x){
      dg(x)
    }
  } else {
    df = function(x){
      dg(x)/g(x,...)
    }
  }
  
  ## Generate starting values  
  x.start = gen.start(x.start,domain)
    
  # Initialization
  set.seed(200)
  x.out=rep(NA,n)
  x.vec = sort(x.start)
  f.vec = f(x.vec)
  if(is.function(dg)){d.vec=df(x.vec)} else {d.vec = grad(f,x.vec)}
  x.temp.start = x.vec[1]
  x.temp.end = x.vec[length(x.vec)]
  dx = (x.temp.end - x.temp.start)/(10*n)
  x.temp = seq(x.temp.start+dx,x.temp.end-dx,by=dx)
  i = 1
  
  # Create vector of samples
  while(i<=n){
    u = runif(1,0,1)
    x.cand = sample(x.temp,size=1, prob=sapply(x.temp, function(x){ 
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
      
      is.bounded(x.cand,x.vec,f.vec,d.vec,f)
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
