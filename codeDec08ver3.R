###################
## Going further ##
###################

# 1. Option to ars() could be log.transformed.f = F
# 2. Derivative of the function can be given by users as: 
#     fprima<-function(x,mu=0,sigma=1){-1/sigma^2*(x-mu)}
#     with parameters by users (?)
# 3. Figure out what to do when x.cand is element in x.vec
# 4. Give the option to give function as input in another way
# 5. Find a way to do the convexity check of f by evaluating fewer f
# 6. Take unbounded domains into consideration
# 7. 

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
# 5. 

################
## User input ##
################

f <- function(x) (-(1/2)*(x)^2)
dom = c(-10,10)
n = 100

###################
## Sub-functions ##
###################

# Derivative of f
df <- function(x) {}
body(df) <- D(body(f), 'x')

# f convexity check
is.convex = function(f,dom){
  x.seq <- seq(dom[1],dom[2], l=100)
  res=T
  for (i in 3:100){
    if((f(x.seq[i])-f(x.seq[i-1])) >= (f(x.seq[i-1])-f(x.seq[i-2]))){
      res = F
    }
  }
  return(res)
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

###################
## Main function ##
###################

ars = function(f, dom, n){
  x.out = rep(NA,n)
  x.vec = c(dom[1],mean(c(dom[1],dom[2])),dom[2])
  f.vec = f(x.vec)
  d.vec = df(x.vec)
  dx = 1e-2
  i = 1
  
  while(i<=n){
    u = runif(1,0,1)
    x.temp = seq(dom[1]+dx,dom[2]-dx,2*dx)
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
      d.vec = append(d.vec,df(x.cand),after=l)
    }
  }
  
  ############
  # test
  u = sapply(x.temp, function(x) upper(x,x.vec,f.vec,d.vec))
  l = sapply(x.temp, function(x) lower(x,x.vec,f.vec))
  plot(c(x.temp,x.temp,x.temp),c(l,u,f(x.temp)))
  ############
  
  return(x.out)
}


##########
## test ##
##########

# User input
f <- function(x) (-(1/2)*(x)^2)
dom = c(-10,10)
n = 100

# Plot f, lower bound and upper bound
dx = 1e-2
x.temp = seq(dom[1]+dx,dom[2]-dx,2*dx)
x.vec = c(dom[1],mean(c(dom[1],dom[2])),dom[2])
f.vec = f(x.vec)
d.vec = df(x.vec)
u = sapply(x.temp, function(x) upper(x,x.vec,f.vec,d.vec))
l = sapply(x.temp, function(x) lower(x,x.vec,f.vec))
plot(c(x.temp,x.temp,x.temp),c(l,u,f(x.temp)))

# Draw a sample using ars()
test.sample = ars(f,dom,n)
hist(test.sample)
