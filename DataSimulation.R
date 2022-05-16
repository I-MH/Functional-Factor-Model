#R code to simulate model (3) on the main paper

require(fda)

FDF.sim <- function(m=50,rangeval=c(0,1),typew=c("MB","PB","Cov"),fn=NULL,epsw=1){
  # it simulates data using model (3) on the main paper
  # args:
  #   m: number of time points to be simulated for a curve 
  #   rangeval: a numeric vector of length 2 defining the interval over which the 
  #              functional data object can be evaluated
  #   typew: type of functional white noise:
  #          MB: Brownian motion
  #          PB: Brownian bridge
  #          Cov: Gaussina process with  covariance Csigma
  #   fn: a matrix with dimension N x K containing the time series \beta_n,k 
  #   epsw: a parameter to modify variance of functional white noise (1=no modification).
  # values: a list  
  #   Data: a matrix of mxN containing the functional data simulated
  #   loadfMat:  a matrix of mxK containing the factor loading curves
  #   tt: points at where the functional data simulated is evaluated
  #   fn: the argument fn  
  
  TT=rangeval[2]
  N <- dim(fn)[1]
  num.f <- dim(fn)[2] # K on the main paper
  typew=match.arg(typew)
  h <- TT/m
  tt<-seq(0,TT,h)
  wn <- fd.wn.sim( N=N,m=m,type = typew,TT=TT )
  base1= create.fourier.basis(nbasis = 5, dropind = 1)
  loadfMat <- eval.basis(evalarg = tt, basisobj = base1 )[,1:num.f] # m x nfactors
  loadfMat <- as.matrix(loadfMat)
  # W_n
  new.wn <-wn
  Datos.f <- loadfMat%*%t(fn) +  epsw*new.wn
  return( list(Data=Datos.f, loadfMat=loadfMat ,tt=tt, fn=fn) )
}


brow<-function(TT,h,N, sigma=1){ 
  # it simulates a brownian motion
  # args:
  #   TT: this is rangeval[2], i.e., (0, TT) is where data is simulated
  #   h: distance of x values where data is sumalated 
  #   N: number of processes to be simulated
  #   sigma: variance parameter. If sigma=.25, the norm of the covariance operator is approx .1
  # values:
  #   Bt: matrix with the values of brownian motions 
  
  m<-length(seq(0,TT,h))
  Aux= matrix(rnorm(n=N*(m-1) , sd=sqrt(h*sigma) ), m-1,N )
  Bt=rbind(rep(0,N),apply(Aux, 2, cumsum)  )
  return(Bt) }


p.brow <- function(TT,h,N,sigma=1){
  # it simulates a brownian bridge
  # args:
  #   TT: this is rangeval[2], i.e., (0, TT) is where data is simulated
  #   h: distance of x values where data is sumalated 
  #   N: number of processes to be simulated
  #   sigma: variance parameter. If sigma=.25, the norm of the covariance operator is approx .1
  # values:
  #   PB: matrix with the values of brownian bridge 
  tt<-seq(0,TT,h)
  n <- length(tt)
  tn <- matrix(tt/tt[n], nrow = n, ncol = N, byrow = FALSE)
  Bt <- brow(TT,h,N,sigma=sigma)
  Aux=matrix( Bt[n,],nrow = n, ncol = N, byrow = TRUE)
  PB <- Bt - tn * Aux
  return(PB) }

CSigma=function(x,y,alpha=.2, beta=.3,norm.op=0.1){
  # it computes the kernel cov, sigma(x,y), values at x, y 
  # this function is used in fd.wn.sim()
  # args:
  #   x: a vector of x values 
  #   y: a vector of y values
  #   alpha, beta: parametres of the cov
  # norm.op: the norm of the cov operator
  # values:
  #   a matrix with the values of sigma(x,y),
  C <- norm.op/0.1816856
  return( outer(x,y, function(x,y) alpha* exp( -beta * abs( x-y ) ) ) )
}

fwn.cov.sim=function(N=100, Cov){
  # it simulates N independent trajectories of Gaussian process with cov=Cov
  # args:
  #   N: number of trajectories 
  #   Cov: a cov matrix
  # values:
  #   a matrix of the values simulated
  nt = ncol(Cov)
  CholCov = chol(Cov)
  return(  t(matrix(rnorm(N * nt), nrow = N, ncol = nt) %*% CholCov) )
}

fd.wn.sim <- function(N=50,m=50,type=c("MB","PB", 'Cov'),TT=1, sigma=1){
  # simulates N functional white noise
  # args:   
  #   N: number of functional white noise
  #   m: number of points to be evaluated for each function
  #   type: type of functional white noise to be simulated
  #     MB: brawnian motion
  #     PB: brawnian bridge
  #     Cov: Gaussina process with  covariance Csigma
  # values:
  #   a matrix of mXN with the values of functions simulated
  dato <- NULL
  h <- TT/m
  tt<-seq(0,TT,h)
  type=match.arg(type)
  if(type=="MB"){
      dato <- brow(TT,h,N, sigma=sigma)  }
  if(type=="PB"){
      dato <- p.brow(TT,h,N) }
  if(type=="Cov"){
    CC=CSigma(x=tt,y=tt)
    dato=white.sim(N=N, Cov=CC)
  }
  return(dato) }

