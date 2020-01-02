## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ------------------------------------------------------------------------
library(Rcpp)
library(GeneralizedHyperbolic)
library(microbenchmark)
source('rwMetropolosR.r')
sourceCpp('rwMetropolosC.cpp')
set.seed(1)
N<-2000
sigma<-c(.05,.5,2,16)
x0<-25
#rwMetropolosR() returns a list with x and k, with the sample vector and the rejection number. 
rwR1<-rwMetropolosR(sigma[1],x0,N)
rwR2<-rwMetropolosR(sigma[2],x0,N)
rwR3<-rwMetropolosR(sigma[3],x0,N)
rwR4<-rwMetropolosR(sigma[4],x0,N)
#rwMetropolosC() returns a vector x, with the first N the sample and the last the rejection number.
rwC1<-rwMetropolosC(sigma[1],x0,N)
rwC2<-rwMetropolosC(sigma[2],x0,N)
rwC3<-rwMetropolosC(sigma[3],x0,N)
rwC4<-rwMetropolosC(sigma[4],x0,N)
refline <- qskewlap(c(.025, .975))
rwR<-list(rwR1,rwR2,rwR3,rwR4)
rwC<-list(rwC1,rwC2,rwC3,rwC4)
##par(mfrow = c(2, 2))
for (i in 1:4){
  plot(1:2000,rwR[[i]]$x,type = "l",xlab = paste("sigma=",sigma[i]),ylab = "x",main = "Using R codes")
  abline(h=refline)
}
for (i in 1:4) {
  plot(1:2000,rwC[[i]][1:2000],type = "l",xlab = paste("sigma=",sigma[i]),ylab = "x",main = "Using Rcpp")
  abline(h=refline)
}
#acceptance rate with rwMetropolosR
print(c(1-rwR1$k/N,1-rwR2$k/N,1-rwR3$k/N,1-rwR4$k/N))
#acceptance rate with rwMetropolosC
print(c(1-rwC1[N+1]/N,1-rwC2[N+1]/N,1-rwC3[N+1]/N,1-rwC4[N+1]/N))


## ------------------------------------------------------------------------
#qqplot
a<-ppoints(100)
##par(mfrow = c(2, 2))
for (i in 1:4) {
  QR<-qskewlap(a)
  Q<-quantile(rwR[[i]]$x,a)
  qqplot(QR,Q,main=paste("Using R codes with sigma=",sigma[i]),xlab = "standard Laplace Quantiles",ylab = "Sample Quantiles")
  lines(QR,QR)
}
for (i in 1:4) {
  QR<-qskewlap(a)
  Q<-quantile(rwC[[i]][1:2000],a)
  qqplot(QR,Q,main=paste("Using Rcpp with sigma=",sigma[i]),,xlab = "standard Laplace Quantiles",ylab = "Sample Quantiles")
  lines(QR,QR)
}


## ------------------------------------------------------------------------
summ<-list(4)
for (i in 1:4) {
  ts<-microbenchmark(rwR[[i]]<-rwMetropolosR(sigma[i],x0,N),rwC[[i]]<-rwMetropolosC(sigma[i],x0,N))
  summ[[i]]<-list(paste0('sigma=',sigma[i]),summary(ts)[,c(1,3,5,6)])
}
print(summ)

