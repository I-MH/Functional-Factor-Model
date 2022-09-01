
# this function requires 
# source('SmoothData.R')
# source('hatK.R')

FDF <- function(argvals=NULL,
                    data, 
                    stationary=TRUE,
                    h=5,
                    k=NULL,
                    kmax=6,
                    p=5,
                    kern_type = "BT",
                    rangeval =c(0,1),
                    nbasis=25,
                    basis='Fourier',
                    lambda = NULL,
                    plot=FALSE){
  # it computes the estimators described on the main paper
  # args:
  #   argvals: a vector containing the points where function values are observed.
  #            If argvals=NULL, this is assumed to be equally spaced in rangeval
  #   data: a matrix with dim m x N, where N is the number of curves
  #   stationary: if TRUE, then it is assumed that data is stationary
  #   h: number of lags to be used
  #   k: number of factor to be used. If NULL, then it is estimated as described on the main paper
  #   kmax: max number of factors to be tested
  #   p: number of eigenfunctions to be used to approximate the inverse of the cov operator
  #   kern_type: type of kernel to be used when estimating the long run cov operator.
  #   rangeval: a numeric vector of length 2 defining the interval over which the functional data object can be evaluated
  #   nbasis: number of basis functions to be used for the functional data.
  #   basis: type of basis. The options are: 'Fourier' and 'Bspline'
  #   lambda: a nonnegative real number specifying the amount of smoothing to be applied 
  #           to the estimated functional parameter. If NULL, this is estimated using
  #            generalized cross validation
  #   plot: if TRUE a plot is displayed
  # values: a list  
  #   hat.beta: the estimated time series
  #   hat.F: the estimated F
  #   hat.K_ratio: the estimated number of factors
  #   Xhat: fitted values of the data
  #   eigenval: eigenvalues of the long run cov 
  #   eigenvalC0: eigenvalues of the cov at lag 0
  #   Ob: (not for users)
  #   Ob_result: (not for users)
  
  kerneltype = switch(kern_type, BT = "Bartlett", PR = "Parzen", 
                      FT = "flat_top", SP = "Simple", flat="flat")
  Kernel <- function(i, h) {
    x = i/h
    if (kerneltype == "flat") {
      return(1)
    }
    if (kerneltype == "Simple") {
      return(0)
    }
    if (kerneltype == "Bartlett") {
      return(1 - x)
    }
    if (kerneltype == "flat_top") {
      if (x < 0.1) {
        return(1)
      }
      else {
        if (x >= 0.1 & x < 1.1) {
          return(1.1 - x)
        }
        else {
          return(0)
        }
      }
    }
    if (kerneltype == "Parzen") {
      if (x < 1/2) {
        return(1 - 6 * x^2 + 6 * abs(x)^3)
      }
      else {
        return(2 * (1 - abs(x))^3)
      }
    }
  }
  
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
  
  CC0 <- t(coef(datafd0))    # N x nbasis
  xbasis <- datafd0$basis
  
  #Stationary case
  if(stationary){
    datafd <- datafd0
    # we need BB, CCmh,CCph
    #BB
    Datopca <- pca.fd(datafd, nharm =p)
    lam0 <- Datopca$values
    bb <- Datopca$harmonics$coefs   # nbasis x num of factors
    BB <-  (1/lam0[1])*bb[,1]%*%t(bb[,1])
    for (l in 2:p) {
      BB <- BB + (1/lam0[l])*bb[,l]%*%t(bb[,l])  }
    # long run cov, this is similar as in LongRun command in fChange package
    CCh <- matrix(0,nbasis,nbasis)
    for (j in 1:h) {
      AuxCCh <- datafd$coefs[,(1+j):N] %*%  t(datafd$coefs[,1:(N-j)]) + 
        t(datafd$coefs[,(1+j):N] %*%  t(datafd$coefs[,1:(N-j)]) )
      CCh <- CCh + Kernel(j, h)*AuxCCh
    }
    # 
    CCh <- (CCh + t(CCh))/2
    #
    # coeff
    CC <- t( coef(datafd) )    # N x nbasis
    #
    Jinprod = inprod(xbasis, xbasis)
    Jsvd <- svd(Jinprod)
    Jaux <- eigen(Jinprod)
    #Jhalf <- Jsvd$u%*%diag(sqrt(Jsvd$d))%*%t(Jsvd$v)
    Jhalf <- Jaux$vectors%*%diag(sqrt(Jaux$values))%*%t(Jaux$vectors)
    Jihalf <- solve(Jhalf)
    # we need eigenvectors of Ob
    Ob <- (1/N)*Jhalf%*%CCh%*%Jinprod%*%BB%*%Jhalf
    rankOb <- qr(Ob)$rank
    result <- eigen(Ob)
    hatK_aux <- hatK(Re(result$values[1:kmax]))
    # num of factors
    if(is.null(k)) k <- hatK_aux$hat_k
    #
    newr <- min(k,rankOb)
    AA <- Jihalf%*%Re(result$vectors[,1:newr])
    result.val <- Re(result$values[1:newr]) #eigenvalues of the operator
    
    #loading factors
    lam <- fd(Re(AA),basisobj = xbasis)
    #factor processes
    ff <- CC%*%Jinprod%*%AA
  }else{
    if(is.numeric(data)){
      newdata <- t(apply(data, 1, diff))
      N <- dim(newdata)[2]
      m <- dim(data)[1]
      datafd <- SmoothData(argvals=argvals,
                           data = newdata, 
                           type_basis=basis, 
                           nbasis=nbasis,
                           rangeval = rangeval,
                           lambda = lambda)$fd
    }else{
      if(sum(class(data)%in%"fd")!=0){
        rangeval=data$basis$rangeval
        newCoeff <- t(apply(data$coefs, 1, diff))
        datafd <- fd(newCoeff, data$basis)
        N <- dim(newCoeff)[2]
        m <- 100
        nbasis <- data$basis$nbasis
      }else{
        stop("data should be a matrix or fd object") 
      }
    }
    
  #BB
  Datopca <- pca.fd(datafd, nharm =p)
  lam0 <- Datopca$values
  bb <- Datopca$harmonics$coefs   # nbasis x num of factors
  BB <-  (1/lam0[1])*bb[,1]%*%t(bb[,1])
  for (l in 2:p) {
    BB <- BB + (1/lam0[l])*bb[,l]%*%t(bb[,l])  }
  # long run cov
  CCh <- matrix(0,nbasis,nbasis)
  for (j in 1:h) {
    AuxCCh <- datafd$coefs[,(1+j):N] %*%  t(datafd$coefs[,1:(N-j)]) + 
      t(datafd$coefs[,(1+j):N] %*%  t(datafd$coefs[,1:(N-j)]) )
    CCh <- CCh + Kernel(j, h)*AuxCCh
  }
  # 
  CCh <- (CCh + t(CCh))/2
  #
  # coeff
  CC <- t( coef(datafd) )    # N x nbasis
  #
  Jinprod = inprod(xbasis, xbasis)
  Jsvd <- svd(Jinprod)
  Jaux <- eigen(Jinprod)
  Jhalf <- Jaux$vectors%*%diag(sqrt(Jaux$values))%*%t(Jaux$vectors)
  Jihalf <- solve(Jhalf)
  #
  Ob <- (1/N)*Jhalf%*%CCh%*%Jinprod%*%BB%*%Jhalf
  rankOb <- qr(Ob)$rank
  result <- eigen(Ob)
  hatK_aux <- hatK(Re(result$values[1:kmax]))
  # num of factors
  if(is.null(k)) k <- hatK_aux$hat_k
  #
  newr <- min(k,rankOb)
  AA <- Jihalf%*%Re(result$vectors[,1:newr])
  result.val <- Re(result$values[1:newr]) 
  
  #loading factors
  lam <- fd(Re(AA),basisobj = xbasis)
  #factor process
  ff <- CC0%*%Jinprod%*%AA
  }
  
  # fitting
  if(is.null(argvals)){
    tt <- seq(rangeval[1],rangeval[2], length.out = m)  
  }else{
    tt=argvals}
  
  Xhat <-  eval.fd(tt,lam)%*%t(as.matrix(ff))
  # scree plot
  hatk0=which.min(Re(result$values[1:kmax]))
  if(plot){
  par(mfrow=c(2,1))
  plot(Re(result$values[1:kmax]),type='o',main = 'Scree Plot', ylab = 'eigenvalues', xlab = 'component number')
  points(hatk0, Re(result$values[1:kmax])[hatk0], col=2,pch=20)
  plot(hatK_aux$ratio,type='o', main = 'Ratio Scree Plot', 
       ylab ='ratio', xlab = 'component number')
  points(hatK_aux$hat_k,hatK_aux$ratio[hatK_aux$hat_k], col=2,pch=20)
  par(mfrow=c(1,1))
  }
  # lam.fd has class pca.fd
  lam.fd <- list(lam, result$values, inprod(datafd, lam))
  class(lam.fd) <- "pca.fd"
  names(lam.fd) <- c("harmonics", "values", "scores")
  ff.Varmx <- varmx.pca.fd(lam.fd, nharm=dim(ff)[2])
  
  return(list(hat.beta=ff, hat.F=lam, Varmx=ff.Varmx,
              hat.beta.Varmx=ff.Varmx$scores, hat.f.Varmx= ff.Varmx$harmonics,
              hat.K_ratio=hatK_aux$hat_k, hat.K.scree=hatk0,
              Xhat=Xhat, eigenval=Re(result$values), eigenvalC0=lam0, Ob=Ob,
              Ob_result=result, argvals=tt))  
}


