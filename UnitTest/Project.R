# Check for concavity across three (x,log(f(x)))
is.local.concave = function(f.vals,x.vals){
  if(((f.vals[3]-f.vals[2])/(x.vals[3]-x.vals[2]))>
     ((f.vals[2]-f.vals[1])/(x.vals[2]-x.vals[1]))){
    stop("Error: The distribution is not log concave!")
  }
  else {
    return("TRUE")
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
