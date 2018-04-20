# angular analysis functions

# A function that calculates the histogram of relative frequencies of points. 

calcHist<-function(x,n=8) {
  h1<-hist((x-(2*pi/(2*n))) %% (2*pi),breaks=seq(0,2*pi,by=2*pi/n),plot=F)$counts
  h1<-h1/sum(h1)
  return(h1)
}

# similar function, but it also considers the order of histogram. Useful for calculating EMD differences.

calcHist1<-function(x,n=8) {
  h1<-hist((x-(2*pi/(2*n))) %% (2*pi),breaks=seq(0,2*pi,by=2*pi/n),plot=F)$counts
  h2<-hist((-x+(2*pi/(2*n))) %% (2*pi),breaks=seq(0,2*pi,by=2*pi/n),plot=F)$counts
  h1<-h1/sum(h1)
  h2<-h2/sum(h2)
  
  if(h2[1]>h1[1]) h1<-h2
  if(h1[8]>h1[2]) {
    h1<-c(h1[1],rev(h1[-1]))
  }
  
  return(h1)
}


# A function that calculates the difference between an input histogram and a histogram of equal frequencies.
calcDiff<-function(x,n=8) {
  a<-rep(1/n,n)
  dif<-sum(abs(x-a))
  return(dif)
}