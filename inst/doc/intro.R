## ------------------------------------------------------------------------
akde<-function(x,h,kernel="gaussian"){
  if (kernel=="gaussian")
    K<-function(u)  {(2*pi)^-0.5*exp(-0.5*u^2)}   else if (kernel=="epanechnikov")
    K<-function(u)  {0.75*(1-u^2)*(abs(u)<=1)}   else if (kernel=="triangular")
    K<-function(u)  {(1-abs(u))*(abs(u)<=1)}    else if (kernel=="triweight")
    K<-function(u)  {35/32*(1-u^2)^3*(abs(u)<=1)}   else if (kernel=="biweight")
    K<-function(u)  {15/16*(1-u^2)^2*(abs(u)<=1)}   else if (kernel=="cosine")
    K<-function(u)  {pi/4*cos(pi*u/2)*(abs(u)<=1)}  else 
    return("ERROR in kernel type");
  n<-length(x)
  ah<-h(x)
  fhat<-function(u){
    sum<-0
    for (i in 1:n) {
      sum<-sum+K((u-x[i])/h(u))/h(u)
    }
    return(sum/n)
  }
  return(list(x=x,kernel=kernel,hx=h,h=ah,fhat=fhat))
}

## ------------------------------------------------------------------------
set.seed(0)
x<-rnorm(2000)
h<-function(u){(10+u^4)^(-0.6)}
fhat<-akde(x,h,kernel = "gaussian")
plot(x,fhat$fhat(x),cex=0.5)

## ------------------------------------------------------------------------
library(FNN)
est.knn<-function(x,k){
  n<-length(x)
  x0<-min(x)
  x1<-max(x)
  f_hat<-function(u){
    kx<-knnx.dist(x,u,k=k)[,k]
    return(k/(2*n*kx)*(u<=x1)*(u>=x0))
  }
  return(f_hat)
}

## ------------------------------------------------------------------------
set.seed(123)
x<-rnorm(1000)
fhat<-est.knn(x,32)
t<-seq(min(x),max(x),length=200)
plot(t,fhat(t),type="l")

## ------------------------------------------------------------------------
set.seed(0)
x<-faithful$eruptions
#h<-function(u){knnx.dist(x,u,k=10)[,k]}
fhat<-akde(x,h=function(u){knnx.dist(x,u,k=10)[,10]},kernel = "gaussian")
t<-seq(min(x),max(x),length=100)
plot(t,fhat$fhat(t),type = "l")

