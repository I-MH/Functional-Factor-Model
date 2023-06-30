# this function requires 
# source('SmoothData.R')
data.fd <- function(argvals=NULL,
                     data,
                     basis,
                     nbasis=25,
                     rangeval =c(0,1),
                     lambda = NULL){
  # returns data with class fd
  if(is.numeric(data)){
    N <- dim(data)[2]
    m <- dim(data)[1]
    datafd0 <- SmoothData(argvals=argvals,
                          data = data, 
                          type_basis=basis, 
                          nbasis=nbasis,
                          rangeval = rangeval,
                          lambda = lambda)$fd 
  }else{ 
    if( sum(class(data)%in%"fd")!=0 ){
      rangeval=data$basis$rangeval
      N <- dim(data$coefs)[2]
      datafd0 <- data
      m <- 100
      nbasis <- data$basis$nbasis
    }else{
      stop("data should be a matrix or fd object") 
    }
  }
  return(list(datafd=datafd0, rangeval=rangeval, N=N,
              m=m, nbasis=nbasis) )
}

diff.fd <- function(datafd){
  # it computes X_n - X_{n-1}
  # datafd= data of class fd from data.fd
  newCoeff <- t(apply(datafd$coefs, 1, diff))
  new.fd <- fd(newCoeff, datafd$basis)
  new.N <- dim(newCoeff)[2]
  datafd$coefs <- newCoeff
  return(datafd)
}


