#' Adaptive Rejection Sampling
#' 
#' @param n Number of samples.
#' @param g Density.
#' @return The result of ars of \code{n} and \code{g}.
#' @examples
#' mysample = ars(100,dnorm,x.start=c(-4,4),plot.type="norm")
#' mysample = ars(100,dnorm,x.start=c(-4,4),plot.type="bounds")
#' mysample = ars(100,dnorm,x.start=c(-4,4),plot.type="acceptance")
#' mysample = ars(100,dnorm,plot.type="bounds")
#' mysample = ars(100,dnorm,plot.type="acceptance")
#' 
#' mysample = ars(n = 100,dbeta, shape1=1.3, shape2=2.7, x.start=c(0.001,0.999), plot.type="bounds")
#' mysample = ars(n = 100,dbeta, shape1=1.3, shape2=2.7, x.start=c(0.001,0.999), plot.type="acceptance")
#' mysample = ars(n = 100,dbeta, shape1=1.3, shape2=2.7, plot.type="bounds", domain=c(0,1))
#' mysample = ars(n = 100,dbeta, shape1=1.3, shape2=2.7, plot.type="acceptance", domain=c(0,1))
#' 
#' f = function(x) {-0.5*x^2}
#' fprima = function(x) {-x}
#' mysample = ars(n=100,f,fprima,log.transform=TRUE,x.start=c(-4,4))

ars = function(n, f, ..., fprima=NULL, log.transform=FALSE, x.start=NULL, 
               domain=c(-Inf,Inf), plot.type = "none",
               pdf.name=""){ 
  
  set.seed(100)
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
    boundary.plot(x.vec,h.vec,d.vec,h=FALSE,pdf.name)
  }
  if(plot.type == "acceptance"){
    boundary.plot(x.vec,h.vec,d.vec,h,pdf.name)
  }
  
  return(x.out)
}
