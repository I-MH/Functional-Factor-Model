library(fChange)

# check this

LongRun.m <-  function (fdobj, h, kern_type = "BT") 
{
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
          s = s + Kernel(i, h) * a
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


