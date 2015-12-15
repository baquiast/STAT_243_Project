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
# +11. Have two starting points
# +12. Vectorize as much as possible
# +13. Efficiency problem

##############
## Comments ##
##############

# 01. "h" is the log of the users density function
# 02. "dom" is the domain defined for the density
# 03. "n" is the requested number of sample elements
# 04. The candidates, x.cand, for a sample value are drawn from the density
#     defined by the upper bound of h. First x.cand is drawn from a equidistant
#     set in dom, with steps 2*dx. The candidate is then adjusted by U(-dx,dx),
#     in order to cover dom.
# 05. Setting plot.type = "acceptance" depend on evaluating h a lot of times and
#     should not be done when h is expensive to evaluate.
# 06. x.start is a vector of starting points
# 07. 

###################
## Sub-functions ##
###################

# Check for log concavity using d(x)
is.local.concave = function(d.vec){
  if(!all.equal(d.vec,sort(d.vec,decreasing=TRUE))){
    stop("Error: The distribution is not log concave!")
  }
}

# Check for log concavity using upper bound, h(x), lower bound 
#is.bounded = function(x.cand, x.vec, h.vec, d.vec, h.cand){
#  u.vec = upper(x.cand, x.vec, h.vec, d.vec)
#  l.vec = lower(x.cand, x.vec, h.vec)
#  if((sum(pmin(u.vec, h.cand, na.rm=FALSE) == h.cand) != length(x.cand))
#     | (sum(pmin(h.cand, l.vec, na.rm=FALSE) == l.vec) != length(x.cand))){
#    stop("Error: The distribution is not log concave!")
#  }
#}

# Evaluate lower bound at x
lower = function(x, x.vec, h.vec){
  z = sort(c(x, x.vec))
  l = match(x,z) - (1:length(x))
  return(((x.vec[l+1]-x) * h.vec[l] + (x-x.vec[l])*h.vec[l+1])/(x.vec[l+1]-x.vec[l]))
}

# Evaluate upper bound at x
upper = function(x, x.vec, h.vec, d.vec){
  z = sort(c(x,x.vec))
  l = match(x,z) - (1:length(x))
  return(pmin(h.vec[l] + (x-x.vec[l])*d.vec[l],
              h.vec[l+1] + (x-x.vec[l+1])*d.vec[l+1] , na.rm=FALSE))
}

# Generate starting points
gen.start = function(x.start,domain,h,fprima,hprima) {
  # When user gives starting values
  if(length(x.start)!=0) {
    # Check if there is at least two unique starting values
    if(length(unique(x.start))<2){
      stop("At least two unique starting values are required")
    }
    # Check if starting values are unique
    #if(length(unique(x.start))!=length(x.start)){
    #  stop("The starting values are not unique")
    #}
    # Check if starting values are inside the domain
    if((sort(domain)[1] > min(x.start)) | (max(x.start) > sort(domain)[2])){
      stop("The starting values are outside the domain")
    }
    x.start = sort(x.start)
  }else{ 
    # When user does not give starting values
    if(length(domain)!=2 | length(unique(domain))!=2){ 
      # Invalid domain
      stop("Invalid domain")
    }else{
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
      p1 = seq(0.1,0,l=10)
      p2 = seq(10,0,l=10)
      p.vec1 = cumsum(exp(p1)/sum(exp(p1)))
      p.vec2 = cumsum(exp(p2)/sum(exp(p2)))
      p.vec = unique(sort(c(p.vec1,p.vec2,0.5*p.vec1,0.5*p.vec2)))
      p.vec = p.vec[p.vec!=0 & p.vec!=1]
      
      x.start = domain
      # Did not vectorize it because it has only two loops
      for(i in 1:2) {
        x.temp1 = rep(NULL,length(p.vec))
        if(h(x.start[i])==Inf){
          x.temp1 = domain[i]-((-1)^(i))*min(1e-2,(domain[2]-domain[1]))*p.vec
          x.temp1 = x.temp1[abs(h(x.temp1))!=Inf]
          x.start[i] = x.temp1[!is.na(grad(h,x.temp1))][1]
        }
        else if(h(x.start[i])==-Inf){
          x.temp1 = domain[i] - ((-1)^(i))*(domain[2]-domain[1])*p.vec
          x.temp1 = x.temp1[abs(h(x.temp1))!=Inf]
          x.start[i] = x.temp1[!is.na(grad(h,x.temp1))][1]
        }
        else {
          if(!is.na(grad(h,x.start[i]))){
            x.start[i] = x.start[i] 
          }else{
            x.temp1 = domain[i]-((-1)^(i))*min(1e-2,(domain[2]-domain[1]))*p.vec
            x.start[i] = x.temp1[!is.na(grad(h,x.temp1))][1]            
          }
        }
      }
      
      #x.mod = optimize(f,interval=x.start,maximum=TRUE)$maximum
      #if(x.start[1]>x.mod | x.mod>x.start[2]) {
      # stop("Cannot find starting points from the given bounds")
      #}
    }
  }
  return(x.start)
}

### The main functions of adaptive rejection sampling
ars.main = function(x.start,h,fprima,hprima,n) {
  
  ## Initialization step
  x.out = NULL
  x.vec = sort(x.start)
  h.vec = h(x.vec)
  if(is.function(fprima)){d.vec=hprima(x.vec)} else {d.vec = grad(h,x.vec)}
  # Evaluate d(x), the derivative of h(x), of starting points    
  xhd.mat = cbind(x.vec,h.vec,d.vec)
  
  x.temp.start = x.vec[1]
  x.temp.end = x.vec[length(x.vec)]
  dx = (x.temp.end-x.temp.start)/(100000000*n)
  x.temp = seq(x.temp.start+dx, x.temp.end-dx,l=20*n)
  
  while(length(x.out) <= n){
    
    ## Sampling step
    u = runif(n-length(x.out),0,1)
    # Sample values from the uniform(0,1) distribution 
    
    x.cand = sample(x.temp,size=n-length(x.out),
                    prob=exp(upper(x.temp,x.vec,h.vec,d.vec)))
    x.cand = runif(n-length(x.out),x.cand-dx,x.cand+dx)
    x.cand = sort(x.cand)
    # Sample values from the s(x), squeezing function  
    
    ## Sampling step: (1) Squeezing test
     # Evaluate l(x) and u(x)
    x.cand.sq1 = x.cand[u <= exp(lower(x.cand,x.vec,h.vec)-
                                   upper(x.cand,x.vec,h.vec,d.vec))]
     # Accepted x.cand
    x.cand.sq0 = x.cand[!x.cand %in% x.cand.sq1]
    x.cand.sq0 = x.cand.sq0[!is.na(x.cand.sq0)] 
    # Rejected x.cand
    
    ## Sampling step: (2) Rejection test
    u.sq0 = u[which(x.cand %in% x.cand.sq0)]
    x.cand.re1 = x.cand.sq0[u.sq0 <= exp(h(x.cand.sq0) - 
                                           upper(x.cand.sq0,x.vec,h.vec,d.vec))] 
    x.cand.re1 = x.cand.re1[!is.na(x.cand.re1)]
     # Accepted x.cand.sq0
    
    x.out = append(x.out, c(x.cand.sq1, x.cand.re1))
     # Accepted samples
    if(length(x.out)==n) {break}
     # Stop the for loop if the number of samples is n
        
    ## Updating step
    h.cand.sq0 = h(x.cand.sq0)
     # Evaluate h(x) of samples rejected at the squeezing test
    
    if(is.function(fprima)){
      d.cand.sq0 = hprima(x.cand.sq0)
    }else{
      d.cand.sq0 = grad(h, x.cand.sq0)
    } # Evaluate d(x) of samples rejected at the squeezing test
    
    cand.mat = cbind(x.cand.sq0, h.cand.sq0, d.cand.sq0)
    
    xhd.mat = cbind(append(xhd.mat[,1],cand.mat[,1]),
                    append(xhd.mat[,2],cand.mat[,2]),
                    append(xhd.mat[,3],cand.mat[,3]))
    xhd.mat = xhd.mat[order(xhd.mat[,1]),]
    x.vec = xhd.mat[,1]
    h.vec = xhd.mat[,2]
    d.vec = xhd.mat[,3]    
    
    is.local.concave(d.vec)  # Check for log concavity
  }
  return(list(xhd.mat,x.out))
}  

# Plot of the acceptance region and bounds
boundary.plot = function(x.vec,h.vec,d.vec,h=FALSE){
  
  # Create upper bound
  nvec = length(x.vec)
  x.vals = seq(x.vec[1],x.vec[nvec],l=100)
  
  x.mat = cbind(x.vec[1:(nvec-1)],x.vec[2:nvec])
  h.mat = cbind(h.vec[1:(nvec-1)],h.vec[2:nvec])
  d.mat = cbind(d.vec[1:(nvec-1)],d.vec[2:nvec])
  x.center = (h.mat[,2]-h.mat[,1]-x.mat[,2]*d.mat[,2]+
                x.mat[,1]*d.mat[,1])/(d.mat[,1]-d.mat[,2])
  h.center = h.mat[,1] + d.mat[,1]*(x.center-x.mat[,1]) 
  x.corners = c(x.vec[1],x.center,x.vec[nvec])
  h.corners = c(h.vec[1],h.center,h.vec[nvec])
  
  # Plot
  if(is.function(h)){
    plot(h, xlim=c(x.vec[1],x.vec[nvec]), lwd =1, yaxt='n', ylab='', xlab="x",
         main ="Acceptance region and rejection of the log density",
         cex.lab = 1.4)
    polygon(x.corners,h.corners,
            col="firebrick", border="firebrick")
    polygon(c(x.vals,x.vec[nvec],x.vec[1]),c(h(x.vals),min(h.vec),min(h.vec)),
            col="palegreen4",border="palegreen4") # acceptance region
    lines(x.corners,h.corners, lwd=1, col="firebrick")
    lines(x.vec,h.vec, lwd =1, col="grey10") # lower bound
    legend(x="topright", legend=c("Acceptance region", "Rejection region",
                                  "Lower bound"),
           col=c("palegreen4","firebrick","grey10"),
           pch=c(15,15,NA),lty=c(0,0,1), seg.len = 0.5, cex = 0.6)
  } else{
    plot(x.corners,h.corners, type = "l",lty =1, lwd=1, col="red", xlab="x",
         yaxt='n', ylab='', main = "Upper and lower bound of the log density",
         cex.lab = 1.4)
    lines(x.vec,h.vec, lwd =1,col="black") # lower bound
    legend(x="topright",legend=c("Upper bound", "Lower bound"),
           col=c("black","red"), lty=c(1,1), cex = 0.6)
  }
}

###################
## Main function ##
###################

ars = function(n, f, ..., fprima=NULL, log.transform=FALSE, x.start=NULL, 
               domain=c(-Inf,Inf), plot.type = "none"){ 
  
  ## Log transformations 
  if(log.transform==TRUE){
    h = function(x){
      if(sum(is.na(f(x,...)))!=0){
        stop("The input density has negative values or too small values")
      }
      return(f(x,...))
    }  
  } else{
    h = function(x){
      if(sum(is.na(log(f(x,...))))!=0){
        stop("The input density has negative values")
      }
      return(log(f(x,...)))
    }
  }
  
  if(log.transform==TRUE){
    hprima = function(x){
      fprima(x)
    }
  } else {
    hprima = function(x){
      fprima(x)/f(x,...)
    }
  }
  
  ## Generate starting points  
  x.start = gen.start(x.start,domain,h,fprima,hprima)
  
  ## Check local concavity for starting points  
  if(is.function(fprima)){d.start=hprima(x.start)} else {d.start = grad(h,x.start)}
  is.local.concave(d.start)
  
  ## ARS main function  
  #set.seed(200)
  xhd.list = ars.main(x.start,h,fprima,hprima,n)
  
  x.vec = xhd.list[[1]][,1]
  h.vec = xhd.list[[1]][,2]
  d.vec = xhd.list[[1]][,3]
  x.out = xhd.list[[2]]
  
  # Optional plot part
  if(plot.type == "bounds"){
    boundary.plot(x.vec,h.vec,d.vec)
  }
  if(plot.type == "acceptance"){
    boundary.plot(x.vec,h.vec,d.vec,h)
  }
  
  return(x.out)
}
