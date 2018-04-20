marpol<-function(mean=0,sd=1) {
  repeat{
    u = runif(1)*2-1
    v = runif(1)*2-1
    s = u*u + v*v
    if(s < 1 & s != 0) break
  }
  mul = sqrt(-2.0*log(s)/s)
  x = mean+sd*u*mul
  y = mean+sd*v*mul
  return(c(x,y))
}

rmarpol<-function(n=100,mean=0,sd=1) {
  return(t(replicate(n,marpol(mean,sd))))
}
