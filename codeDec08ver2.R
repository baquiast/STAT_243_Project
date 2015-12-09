# If this function is given.
f <- function(x) (-(1/2)*(x)^2) #### If a user gives a function that is not log transformed?

# First, we can get the derivative of the function
g <- function(x) {} # Change the function name!
body(g) <- D(body(f), 'x')
# Or, the derivative of the function can be also given by users as: fprima<-function(x,mu=0,sigma=1){-1/sigma^2*(x-mu)}
 # put parameters by users (?)



# log convexity check
 # Question: Should I write a code to log transform any given function? 
rand.x <- runif(n=100,min=-100,max=100) # Random sample x's from unif distribution ===> Don't have to be random!!!!! 
rand.x <- sort(rand.x)                  # Sort the x's  
rand.fx <- f(rand.x)                    # Obtain the f(x)'s (Note that f(x) is log-transformed distribution)
rand.sl <- rep(NA,99)                   # Make an empty vector to store the values of slopes
for (i in 2:100) {                      # Compute slopes of the function and store them in the vector 
  rand.sl[i-1] <- (rand.fx[i]-rand.fx[i-1])/(rand.x[i]-rand.x[i-1])
}
rand.slincr <- rep(NA,98)               # Make an empty vector to store the values of the changes in the slopes
for (i in 3:100) {                      # Compute the changes in the slopes and store them 
  rand.slincr[i-2] <- rand.sl[i-1]-rand.sl[i-2]
}
if (max(rand.slincr)>=0) {              # Note that distributions are not log convex if the slope changes are nonnegative
  print("This distribution is not log convex")
} 



  
# adaptive rejection sampling
 # 1) Initialization
n <- 100
x.out <- rep(NA,n) # If I want to sample out 100 x's
 # 1-1) Select the three starting points => ########## domain as input?
start.select=0                       
while (start.select==2) {
  init <- runif(n=3,min=-n,max=n) 
  # Question: -100, 100 are reasonable minimum and maximum?
  init <- sort(init)
  start.select <- (f(init[2])-f(init[1])>0)+(f(init[3])-f(init[2])<0)
}
x.out[1:3] <- init

# 1-2) Make upper bounds
u.bound1 <- function(x) (f(x.out[1])+g(x.out[1])*(x-x.out[1])) # from x1 to u.inter12
u.bound2 <- function(x) (f(x.out[2])+g(x.out[2])*(x-x.out[2])) # from u.inter12 to u.inter23
u.bound3 <- function(x) (f(x.out[3])+g(x.out[3])*(x-x.out[3])) # from u.inter23 to x3
u.inter12 <- function(x) u.bound1(x) - u.bound2(x)
u.inter23 <- function(x) u.bound2(x) - u.bound3(x)
u.inter12 <- uniroot(u.inter12,c(x.out[1],x.out[2]))$root
u.inter23 <- uniroot(u.inter23,c(x.out[2],x.out[3]))$root
 # 1-3) Make lower bounds
l.bound1 <- function(x) (f(x.out[1])+((f(x.out[2])-f(x.out[1]))/(x.out[2]-x.out[1]))*(x-x.out[1])) # from x1 to x2
l.bound2 <- function(x) (f(x.out[2])+((f(x.out[3])-f(x.out[2]))/(x.out[3]-x.out[2]))*(x-x.out[2])) # from x2 to x3
 # 1-4) Sampling out an x candidate from the distribution of the upper bound
rand.xden <- cbind(x=seq(x.out[1],x.out[3],length.out=100),den=rep(0,100)) # 100 is reasonable?
for (i in 1:100) {
  if (rand.xden[i,1] < u.inter12) {rand.xden[i,2] <- exp(u.bound1(rand.xden[i,1]))
   } else if (rand.xden[i,1] < u.inter23) {rand.xden[i,2]<- exp(u.bound2(rand.xden[i,1]))
   } else {rand.xden[i,2] <- exp(u.bound3(rand.xden[i,1]))}
} ## right?
x.cand <- sample(rand.xden[,1],size=1,prob=rand.xden[,2])

 # 2) Sampling step
unif.rand <- runif(n=1,min=0,max=1)

if (x.cand < x[2]) {l.max <- l.bound1(x.cand)
 } else {l.max <- l.bound2(x.cand)}

if (x.cand < u.inter12) {u.max <- u.bound1(x.cand)
 } else if (x.cand < u.inter23) {u.max <- u.bound2(x.cand)
 } else {u.max <- u.bound3(x.cand)}

 # 3) Updating step

x.sample <- x.out

if (u <= exp(l.max-u.max)) {
  x.out[4] <- x.cand
  } else {
  d.max <- f(x.cand)
  if (u > exp(d.max-u.max)) {
    x.out[4] <- x.cand
  }
  x.sample[4] <- x.cand
  x.sample <- sort(x.sample) # First sort the x's again to make bounds
  
  u.bound1 <- function(x) (f(x.sample[1])+g(x.sample[1])*(x-x.sample[1])) # from x1 to u.inter12
  u.bound2 <- function(x) (f(x.sample[2])+g(x.sample[2])*(x-x.sample[2])) # from u.inter12 to u.inter23
  u.bound3 <- function(x) (f(x.sample[3])+g(x.sample[3])*(x-x.sample[3])) # from u.inter23 to x3
  u.bound4 <- function(x) (f(x.sample[4])+g(x.sample[4])*(x-x.sample[4])) # from u.inter34 to x4
  u.inter12 <- function(x) u.bound1(x) - u.bound2(x)
  u.inter23 <- function(x) u.bound2(x) - u.bound3(x)
  u.inter34 <- function(x) u.bound3(x) - u.bound4(x)
  u.inter12 <- uniroot(u.inter12,c(x.sample[1],x.sample[2]))$root
  u.inter23 <- uniroot(u.inter23,c(x.sample[2],x.sample[3]))$root
  u.inter34 <- uniroot(u.inter34,c(x.sample[3],x.sample[4]))$root
  
  l.bound1 <- function(x) (f(x.sample[1])+((f(x.sample[2])-f(x.sample[1]))/(x.sample[2]-x.sample[1]))*(x-x.sample[1])) # from x1 to x2
  l.bound2 <- function(x) (f(x.sample[2])+((f(x.sample[3])-f(x.sample[2]))/(x.sample[3]-x.sample[2]))*(x-x.sample[2])) # from x2 to x3
  l.bound3 <- function(x) (f(x.sample[3])+((f(x.sample[4])-f(x.sample[3]))/(x.sample[4]-x.sample[3]))*(x-x.sample[3])) # from x3 to x4
} 