test.ars.main = function() {
  x.start <- c(-4,4)
  h = function(x) {log(dnorm(x,0,1))}
  fprima = NULL
  hprima = NULL
  n=10
  list <- ars.main(x.start,h,fprima,hprima,n)
  checkEquals(10,length(list[[2]]))
}
