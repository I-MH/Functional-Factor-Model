# computes Long run cov (several versions)

Kernel <- function(i, h, kerneltype) {
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


lrc <- function(datafd, kern_type,h){
  # returns the eigenfunctions of the long run cov
  kerneltype = switch(kern_type, BT = "Bartlett", PR = "Parzen", 
                      FT = "flat_top", SP = "Simple", flat="flat")
  
  N <- dim(datafd$coefs)[2]
  nbasis <- datafd$basis$nbasis
  CCh0 <- datafd$coefs[,1:N]%*%t(datafd$coefs[,1:N])
  CCh <- matrix(0,nbasis,nbasis) + CCh0
  for (j in 1:h) {
    AuxCCh <- datafd$coefs[,(1+j):N]%*%t(datafd$coefs[,1:(N-j)])+ 
      t( datafd$coefs[,(1+j):N] %*%  t(datafd$coefs[,1:(N-j)]) )
    CCh <- CCh + Kernel(j, h,kerneltype)*AuxCCh
  }
  CCh <- (CCh + t(CCh))/2 # in case is not symmetric
  
  xbasis <- datafd$basis
  Jinprod = inprod(xbasis, xbasis)
  Jaux <- eigen(Jinprod)
  Jhalf <- Jaux$vectors%*%diag(sqrt(Jaux$values))%*%t(Jaux$vectors)
  Jihalf <- solve(Jhalf)
  Ob <- (1/N)*Jhalf%*%CCh%*%Jhalf
  result <- eigen(Ob)
  bb <- Jihalf%*%Re(result$vectors) # coeff of the eigenfunctions
  eigenval <- Re(result$values) #eigenvalues of the operator
  eigenf <- fd(Re(bb), basisobj = xbasis)
  return(list(lrc.k=CCh, Ob=Ob, bb=Re(bb), 
              eigenf=eigenf, eigenval=eigenval))
}

lrc.paper <- function(datafd, kern_type,h,p){
  # returns the eigenfunctions of the long run cov
  kerneltype = switch(kern_type, BT = "Bartlett", PR = "Parzen", 
                      FT = "flat_top", SP = "Simple", flat="flat")
  
  N <- dim(datafd$coefs)[2]
  nbasis <- datafd$basis$nbasis
  
  Datopca <- pca.fd(datafd, nharm =p)
  lam0 <- Datopca$values
  bb0 <- Datopca$harmonics$coefs   # nbasis x num of factors
  BB <-  (1/lam0[1])*bb0[,1]%*%t(bb0[,1])
  for (l in 2:p) {
    BB <- BB + (1/lam0[l])*bb0[,l]%*%t(bb0[,l])  }
  # long run cov, this is similar as in LongRun command in fChange package
  CCh <- matrix(0,nbasis,nbasis)
  for (j in 1:h) {
    AuxCCh <- datafd$coefs[,(1+j):N] %*%  t(datafd$coefs[,1:(N-j)]) + 
      t(datafd$coefs[,(1+j):N] %*%  t(datafd$coefs[,1:(N-j)]) )
    CCh <- CCh + Kernel(j, h, kerneltype)*AuxCCh
  }
  # 
  CCh <- (CCh + t(CCh))/2 # in case is not symmetric
  
  xbasis <- datafd$basis
  Jinprod = inprod(xbasis, xbasis)
  Jaux <- eigen(Jinprod)
  Jhalf <- Jaux$vectors%*%diag(sqrt(Jaux$values))%*%t(Jaux$vectors)
  Jihalf <- solve(Jhalf)
  Ob <- (1/N)*Jhalf%*%CCh%*%Jinprod%*%BB%*%Jhalf
  result <- eigen(Ob)
  bb <- Jihalf%*%Re(result$vectors) # coeff of the eigenfunctions
  eigenval <- Re(result$values) #eigenvalues of the operator
  eigenf <- fd(Re(bb), basisobj = xbasis)
  return(list(lrc.k=CCh, Ob=Ob, bb=Re(bb), 
              eigenf=eigenf, eigenval=eigenval))
}

lrc.paper.NoInv <- function(datafd, kern_type,h, p){
  # returns the eigenfunctions of the long run cov
  kerneltype = switch(kern_type, BT = "Bartlett", PR = "Parzen", 
                      FT = "flat_top", SP = "Simple", flat="flat")
  
  N <- dim(datafd$coefs)[2]
  nbasis <- datafd$basis$nbasis
  CCh <- matrix(0,nbasis,nbasis)
  for (j in 1:h) {
    AuxCCh <- datafd$coefs[,(1+j):N]%*%t(datafd$coefs[,1:(N-j)])+ 
      t( datafd$coefs[,(1+j):N] %*%  t(datafd$coefs[,1:(N-j)]) )
    CCh <- CCh + Kernel(j, h,kerneltype)*AuxCCh
  }
  CCh <- (CCh + t(CCh))/2 # in case is not symmetric
  
  xbasis <- datafd$basis
  Jinprod = inprod(xbasis, xbasis)
  Jaux <- eigen(Jinprod)
  Jhalf <- Jaux$vectors%*%diag(sqrt(Jaux$values))%*%t(Jaux$vectors)
  Jihalf <- solve(Jhalf)
  Ob <- (1/N)*Jhalf%*%CCh%*%Jhalf
  result <- eigen(Ob)
  bb <- Jihalf%*%Re(result$vectors) # coeff of the eigenfunctions
  eigenval <- Re(result$values) #eigenvalues of the operator
  eigenf <- fd(Re(bb), basisobj = xbasis)
  return(list(lrc.k=CCh, Ob=Ob, bb=Re(bb), 
              eigenf=eigenf, eigenval=eigenval))
}

lrc.fchange <- function(datafd, kern_type,h, p){
  # returns the eigenfunctions of the long run cov
  N <- dim(datafd$coefs)[2]
  nbasis <- datafd$basis$nbasis
  
  LR <- LongRun.m(fdobj = datafd, h=h, 
                 kern_type = kern_type)
  CCh <- LR$covm
  xbasis <- datafd$basis
  Jinprod = inprod(xbasis, xbasis)
  Jaux <- eigen(Jinprod)
  Jhalf <- Jaux$vectors%*%diag(sqrt(Jaux$values))%*%t(Jaux$vectors)
  Jihalf <- solve(Jhalf)
  Ob <- (1/N)*Jhalf%*%CCh%*%Jhalf
  result <- eigen(Ob)
  bb <- Jihalf%*%Re(result$vectors) # coeff of the eigenfunctions
  eigenval <- Re(result$values) #eigenvalues of the operator
  eigenf <- fd(Re(bb), basisobj = xbasis)
  return(list(lrc.k=CCh, Ob=Ob, bb=Re(bb), 
              eigenf=eigenf, eigenval=eigenval))
}

LongRun.m <-  function (fdobj, h, kern_type) { 
  # used in lrc.fchange()
  kerneltype = switch(kern_type, BT = "Bartlett", PR = "Parzen", 
                      FT = "flat_top", SP = "Simple", flat="flat")
  
  N = ncol(fdobj$coefs)
  D = nrow(fdobj$coefs)
  basis = fdobj$basis
  D_mat = matrix(0, D, D)
  fdobj_centered = center.fd(fdobj)
  for (k in 1:D) {
    for (r in k:D) {
      s = fdobj_centered$coefs[k, 1:N] %*% fdobj_centered$coefs[r, 
                                                                1:N]
      if (h > 0) {
        for (i in 1:h) {
          a = fdobj_centered$coefs[k, 1:(N - i)] %*% 
            fdobj_centered$coefs[r, (i + 1):N]
          a = a + fdobj_centered$coefs[r, 1:(N - i)] %*% 
            fdobj_centered$coefs[k, (i + 1):N]
          s = s + Kernel(i, h, kerneltype) * a
        }
      }
      D_mat[k, r] = s
      D_mat[r, k] = D_mat[k, r]
    }
  }
  D_mat = D_mat/N
  eigen_struct = eigen(D_mat, symmetric = TRUE)
  eigenfunc = fd(eigen_struct$vectors, basisobj = basis)
  
  list(e_fun = eigenfunc, e_val = abs(eigen_struct$values), 
       covm = D_mat)
  
}
