set.seed(200)

The normal distribution and multivariate normal distributions.
N(0,1)  [-4,4]
F_N = function(x) {dnorm(x,0,1)}
plot(F_N,-4,4)
x_N <- c(-4,0,4)
f_N <- log(F_N(x_N))
is.local.concave(f_N,x_N)
mysample = ars(5,dnorm,x.start=c(-4,0,4),plot.type="bounds")
u      0.5337724
x.cand -1.12
x.cand -1.091335
u      0.6910399
x.cand 1.44
x.cand 1.548574
x.vec -4.000000 -1.091335  0.000000  1.548574  4.000000
f.vec -8.9189385 -1.5144445 -0.9189385 -2.1179792 -8.9189385


The exponential distribution.
Exp(1)  [0,4]
F_Exp = function(x) {dexp(x,1)}
plot(F_Exp,0,4)
x_Exp <- c(0,4,100)
f_Exp <- log(F_Exp(x_Exp))
is.local.concave(f_Exp,x_Exp)


The uniform distribution over any convex set.
U[0,1]  [0,1]
F_U = function(x) {dunif(x,0,1)}
plot(F_U,0,1)
x_U <- c(0,0.5,1)
f_U <- log(F_U(x_U))
is.local.concave(f_U,x_U)

The logistic distribution.
Log[0,1]  [-4,4]
F_Log = function(x) {dlogis(x,0,1)}
plot(F_Log,-4,4)
x_Log <- c(-4,0,4)
f_Log <- log(F_Log(x_Log))
is.local.concave(f_Log,x_Log)

The extreme value distribution.


The Laplace distribution.
Lap[0,1]  [-4,4]
library(VGAM)
F_Lap = function(x) {dlaplace(x,0,1)}
plot(F_Lap,-4,4)
x_Lap <- c(-4,0,4)
f_Lap <- log(F_Lap(x_Lap))
is.local.concave(f_Lap,x_Lap)

The chi distribution.
Chi(2)  [0,10]
F_Chi = function(x) {dchi(x,2)}
plot(F_Chi,0,10)
x_Chi <- c(0,4,10)
f_Chi <- log(F_Chi(x_Chi))
is.local.concave(f_Chi,x_Chi)

The Wishart distribution, where n >= p + 1.[3]

The Dirichlet distribution, where all parameters are >= 1.[3]
Dirichlet(1,2,3)  [0,1]
library(fields)
library(gtools)
library(NCmisc)
F_D = function(x) {ddirichlet(x,c(1,2,3))}
x_D <- matrix(nrow=3,ncol=3)
f_D <- c()
for(i in 1:3) {
x_D[i,] <- c(i/10,2*i/10,1-(3*i/10))
f_D[i] <- log(F_D(x_D[i,]))
}
is.local.concave(f_D,x_D)

The gamma distribution if the shape parameter is >= 1.
Gamma(1,2)  (0,4]
F_G = function(x) {dgamma(x,1,2)}
plot(F_G,0,4)
x_G <- c(1,2,4)
f_G <- log(F_G(x_G))
is.local.concave(f_G,x_G)

The chi-square distribution if the number of degrees of freedom is >= 2.
Chi-Square(2)  [0,10]
F_ChiSq = function(x) {dchisq(x,2)}
plot(F_ChiSq,0,10)
x_ChiSq <- c(0,4,10)
f_ChiSq <- log(F_ChiSq(x_ChiSq))
is.local.concave(f_ChiSq,x_ChiSq)

The beta distribution if both shape parameters are >= 1.
Beta(3,3)  [0,1]
F_B = function(x) {dbeta(x,3,3)}
plot(F_B,0,1)
x_B <- c(0.2,0.4,0.9)
f_B <- log(F_B(x_B))
is.local.concave(f_B,x_B)

The Weibull distribution if the shape parameter is >= 1.
Weibull(2,3)  [0,10]
F_W = function(x) {dweibull(x,2,3)}
plot(F_W,0,10)
x_W <- c(0.5,2,5)
f_W <- log(F_W(x_W))
is.local.concave(f_W,x_W)

The Student t-distribution.
t(2)   [-10,10]
F_t = function(x) {dt(x,2)}
plot(F_t,-10,10)
x_t <- c(0,1,10)
f_t <- log(F_t(x_t))
is.local.concave(f_t,x_t)

The Cauchy distribution.
Cauchy(1,0.5) [0,4]
F_C = function(x) {dcauchy(x,1,0.5)}
plot(F_C,0,4)
x_C <- c(0,0.2,0.8)
f_C <- log(F_C(x_C))
is.local.concave(f_C,x_C)

The Pareto distribution.


The log-normal distribution.
Lognormal(0,1)
F_L = function(x) {dlnorm(x)}
plot(F_L,0,4)
x_L <- c(0.1,1,100)
f_L <- log(F_L(x_L))
is.local.concave(f_L,x_L)

The F-distribution.
F(1,2) [0,4]
F_F = function(x) {df(x,1,2)}
plot(F_F,0,4)
x_F <- c(0,0.2,0.8)
f_F <- log(F_F(x_F))
is.local.concave(f_F,x_F)
