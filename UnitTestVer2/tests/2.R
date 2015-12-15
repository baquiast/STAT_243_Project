test.lower <- function() {
  x.vec <- c(-4,-1.091335,0,1.548574,4)
  f.vec <- c(-8.9189385,-1.5144445,-0.9189385,-2.1179792,-8.9189385)
  checkEquals(lower(1.548574,x.vec,f.vec),-2.117979,tolerance=10^-6)
  checkEquals(lower(-1.091335,x.vec,f.vec),-1.514444,tolerance=10^-6)
  checkEquals(lower(0,x.vec,f.vec),-0.9189385,tolerance=10^-6)
  checkEquals(lower(4,x.vec,f.vec),-8.9189385,tolerance=10^-6)
  checkEquals(lower(-2,x.vec,f.vec),-3.827604,tolerance=10^-6)
  checkEquals(lower(2,x.vec,f.vec),-3.370365,tolerance=10^-6)
}

test.upper <- function() {
  x.vec <- c(-4,-1.091335,0,1.548574,4)
  f.vec <- c(8.9189385,-1.5144445,-0.9189385,-2.1179792,-8.9189385)
  d.vec <- c(4,1.091335,0,-1.548574,-4)
  checkEquals(upper(1.548574,x.vec,f.vec,d.vec),-2.117979,tolerance=10^-6)
  checkEquals(upper(-1.091335,x.vec,f.vec,d.vec),-1.514444,tolerance=10^-6)
  checkEquals(upper(0,x.vec,f.vec,d.vec),-0.918939,tolerance=10^-6)
  checkEquals(upper(4,x.vec,f.vec,d.vec),-8.918939,tolerance=10^-6)
  checkEquals(upper(-2,x.vec,f.vec,d.vec),-2.506102,tolerance=10^-6)
  checkEquals(upper(2,x.vec,f.vec,d.vec),-2.817046,tolerance=10^-6)
}