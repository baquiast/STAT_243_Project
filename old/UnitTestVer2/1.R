is.local.concave = function(d.vec){
  if(!all.equal(d.vec,sort(d.vec,decreasing=TRUE))==TRUE){
    stop("Error: The distribution is not log concave!")
  } else {
    return("TRUE")
  }
}

is.bounded = function(x.cand.re0, x.vec, f.vec, d.vec, f.cand.re0){
  if(sum(pmin(upper(x.cand.re0, x.vec, f.vec, d.vec),f.cand.re0, na.rm=FALSE)==f.cand.re0)!=length(x.cand.re0)
     | sum(pmin(f.cand.re0,lower(x.cand.re0, x.vec, f.vec),na.rm=FALSE)==lower(x.cand.re0, x.vec, f.vec))!=length(x.cand.re0)) {
    stop("Error: The distribution is not bounded!")
  } else {
    return("TRUE")
  }
}

# Evaluate lower bound at x
lower = function(x, x.vec, f.vec){
  z = sort(append(x,x.vec))
  l = match(x,z)-1:length(x)
  return(((x.vec[l+1]-x)*f.vec[l]+(x-x.vec[l])*f.vec[l+1])/(x.vec[l+1]-x.vec[l]))
}

# Evaluate upper bound at x
upper = function(x, x.vec, f.vec, d.vec){
  z = sort(append(x,x.vec))
  l = match(x,z)-1:length(x)
  return(pmin(f.vec[l]+(x-x.vec[l])*d.vec[l],
              f.vec[l+1]+(x-x.vec[l+1])*d.vec[l+1],na.rm=FALSE))
}

# Generate starting points
gen.start = function(x.start,domain,f) {
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
