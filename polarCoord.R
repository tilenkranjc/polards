cart2pol <- function(x,y=0) {
  if(is.vector(x) & length(x)==2 & y == 0) {
    y<-x[2]
    x<-x[1]
  } else if(!is.vector(x) & length(x) > 1) {
    if(nrow(x)>2 & y == 0) {
      y<-x[,2]
      x<-x[,1]
    }
  }
  b <- complex(real=x, imaginary = y)
  data.frame(radius = Mod(b), angle=Arg(b))
}

meanAngle<-function(x) {
  m<-Arg(length(x)*sum(exp(complex(imaginary=x))))
  return(m)
}