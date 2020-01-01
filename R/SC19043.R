#' @title An adaptive KDE
#' @description An adaptive kernel density estimator using R.
#' @param x the data from which the estimate is to be computed.
#' @param h the bandwidth function to be used.
#' @param kernel a character string giving the smoothing kernel to be used, with default "gaussian".
#' @return x  data points - same as input. 
#' @return kernel name of kernel to use.
#' @return hx the bandwidth function - same as input.
#' @return h the bandwidth value at x to use.
#' @return fhat the estimated density function. 
#' @examples
#' \dontrun{
#' x<-rnorm(2000)
#' h<-function(u){(10+u^4)^(-0.6)}
#' fhat<-akde(x,h,kernel="gaussian")
#' plot(x,fhat$fhat(x),cex=0.5)
#' }
#' @export
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


#' @title KNN density estimator
#' @description KNN density estimator
#' @param x the data from which the estimate is to be computed.
#' @param k the maximum number of nearest neighbors to search, usually close to k^-0.5.
#' @return the estimated density function. 
#' @import "FNN"
#' @examples
#' \dontrun{
#' x<-rnorm(1000)
#' fhat<-est.knn(x,32)
#' t<-seq(min(x),max(x),length=200)
#' plot(t,fhat(t),type="l")
#' }
#' @export
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



