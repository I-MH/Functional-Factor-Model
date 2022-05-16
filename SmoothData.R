
require(fda)
SmoothData <- function(argvals=NULL, 
                       data,
                       type_basis=c('Fourier','Bspline'),
                       nbasis=25,
                       rangeval =c(0,1),
                       lambda=NULL){
  # this function smooths data with fda package using a basis function.
  # Lambda is selected such that it minimizes the generalized cross validation
  #
  # args:
  #   argvals: a vector containing the points where function values are observed.
  #            If argvals=NULL, this is assumed to be equally spaced in rangeval
  #   data: a matrix with dim m x N, where N is the number of curves 
  # values: fda object
  
  if(is.null(argvals)){
    tt <- seq(rangeval[1],rangeval[2], length.out = dim(data)[1])  
  }else{
    tt=argvals }
  
  type_basis= match.arg(type_basis)
  
  if(type_basis=='Fourier'){
    base <- create.fourier.basis(rangeval =rangeval, nbasis=nbasis)
  }else{
    base <- create.bspline.basis(rangeval = rangeval, nbasis =  nbasis, norder = 4) }
  
  if(is.null(lambda)){
    lambda=1e-6
    Lfdobj <- 2
    datodPar <- fdPar(base, Lfdobj, lambda = lambda)
    # selecting the best lambda
    lambdas = 10^seq(-8,3,by=0.5)    # possoble options for lambda
    mean.gcv = rep(0,length(lambdas)) # store mean gcv
    for(ilam in 1:length(lambdas)){
      # Set lambda
      curv.fdPari = datodPar
      curv.fdPari$lambda = lambdas[ilam]
      # Smooth
      tempSmoothi = smooth.basis(tt,data,curv.fdPari)
      # Record average gcv
      mean.gcv[ilam] = mean(tempSmoothi$gcv)
    }
    #select the lowest of these and smooth
    best = which.min(mean.gcv)
    lambdabest = lambdas[best]
    datodPar$lambda = lambdabest
    datafd = smooth.basis(tt,data,datodPar)
  }else{
    Lfdobj <- 2
    datodPar <- fdPar(base, Lfdobj, lambda = lambda)
    datafd <- smooth.basis(tt, data, datodPar)
  }
  return(datafd)
}
  
  
  