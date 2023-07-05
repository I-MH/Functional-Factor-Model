
# this function requires 
# source('SmoothData.R')
# source('hatK.R')
# source('data.fd.R')
# source('lrc.R')

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
                 type_lrc=c('paper','paperNoInv', 'lrc', 
                            'lrcFc'),
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
  
  # transform data to fd class ------------------------------
  newData <- data.fd(argvals=argvals,
                     data=data,
                     basis=basis,
                     nbasis=nbasis,
                     rangeval =rangeval,
                     lambda = lambda)
    
  datafd0 <- newData$datafd
  rangeval <- newData$rangeval
  N <- newData$N
  m <- newData$m
  nbasis <- newData$nbasis
  CoeffData <- t(coef(datafd0))    # N x nbasis
  xbasis <- datafd0$basis
  Jinprod = inprod(xbasis, xbasis)

  if(stationary){
    datafd <- datafd0
    }else{
    datafd <- diff.fd(datafd = datafd0)
  }
  
  # estimate long run cov ----------------------------------- 
  type_lrc= match.arg(type_lrc)
  if(type_lrc=='paper'){
    result.lrc <- lrc.paper(datafd=datafd, 
                      kern_type = kern_type, 
                      h=h, p=p)  
  }else{
    if(type_lrc=='paperNoInv'){
      result.lrc <- lrc.paper.NoInv(datafd=datafd, 
                              kern_type = kern_type, 
                              h=h)  
    }else{
      if(type_lrc=='lrc'){
        result.lrc <- lrc(datafd=datafd, 
                          kern_type = kern_type, 
                          h=h) 
      }else{
        # this is lrc obtained from fChange package
        result.lrc <- lrc.fchange(datafd=datafd, 
                          kern_type = kern_type, 
                          h=h, p=p) 
      }
    }
  }
  
  # estimate num of factors ---------------------------------
  hatK_aux <- hatK(Re(result.lrc$eigenval[1:kmax]))
  if(is.null(k)){
    k <- hatK_aux$hat_k
  } 
  rank_bb <- qr(result.lrc$Ob)$rank 
  newr <- min(k,rank_bb)
  # estimating factors ---------------------------
  bb.f <- result.lrc$bb[,1:newr] # coef for factor
  e.val <- Re(result.lrc$eigenval[1:newr]) #eigenvalues of the operator
  #loading factors
  lam <- fd(Re(bb.f), basisobj = xbasis)
  #factor processes
  ff <- CoeffData%*%Jinprod%*%bb.f 
  
  # fitting -------------------------------------
  if(is.null(argvals)){
    tt <- seq(rangeval[1],rangeval[2], length.out = m)  
  }else{
    tt=argvals}
  Xhat <-  eval.fd(tt,lam)%*%t(as.matrix(ff))
  # fd version
  coeff.fit <- lam$coefs %*% t(ff)
  Xhat.fd <- fd(coeff.fit, basisobj = xbasis)
  error.fd <- fd(datafd0$coefs-coeff.fit, basisobj = xbasis)
  
  # obtain scree plot -----------------------------
  hatk0=which.min(Re(result.lrc$eigenval[1:kmax]))
  if(plot){
    par(mfrow=c(2,1))
    plot(Re(result.lrc$eigenval[1:kmax]),
         type='o',main = 'Scree Plot', 
         ylab = 'eigenvalues', xlab = 'component number')
    points(hatk0, 
           Re(result.lrc$eigenval[1:kmax])[hatk0], 
           col=2,pch=20)
    plot(hatK_aux$ratio,type='o', main = 'Ratio Scree Plot', 
         ylab ='ratio', xlab = 'component number')
    points(hatK_aux$hat_k,hatK_aux$ratio[hatK_aux$hat_k], col=2,pch=20)
    par(mfrow=c(1,1))
  }
  
  # make lam.fd a pca.fd class ------------------------------- 
  # this is to for varmx.pca.fd()
  lam.fd <- list(lam, result.lrc$eigenval, ff)
  class(lam.fd) <- "pca.fd"
  names(lam.fd) <- c("harmonics", "values", "scores")
  ff.Varmx <- varmx.pca.fd(lam.fd, nharm=dim(ff)[2])
  
  return(list(hat.beta=ff, hat.F=lam, Varmx=ff.Varmx,
              hat.beta.Varmx=ff.Varmx$scores, hat.f.Varmx= ff.Varmx$harmonics,
              hat.K_ratio=hatK_aux$hat_k, hat.K.scree=hatk0, h=h, type_lrc=type_lrc,
              kern_type=kern_type, Xhat=Xhat, Xhat.fd=Xhat.fd, 
              error.fd=error.fd, eigenval=Re(result.lrc$eigenval), argvals=tt))  
}


