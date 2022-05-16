
source('Est_FDF.R')
source('hatK.R')
source('SmoothData.R')

# laod data FedYieldCurve.RData

timeYear <- c(3/12, 6/12,1,2,3,5,7,10) # time in years
tt <- c( (3/12)/10, (6/12)/10,1/10,2/10,3/10,5/10,7/10,10/10) # in (0,1)
N <- dim(MatYield)[2] # sample size 

# plot data
matplot(timeYear, MatYield, type = 'l', xlab = 'maturities structure', 
        ylab = '', main= "Yield Curves")

# test for stationarity
library(ftsa)
T_stationary(MatYield, L=15) # not stationary

# estimate parameters for non stationary components
est_nst <- Est_FDF(argvals=tt, 
                      data=MatYield,
                      h=2,k=1,p=4,
                      nbasis=15,
                      basis='Bspline',
                      kern_type = "BT",
                      lambda = NULL,
                      stationary = FALSE,
                      plot = FALSE)

# subtract the nonstationary component
MatYield2 <- MatYield - est_nst$Xhat
# test
T_stationary(MatYield2, L=15)

# estimate parameters for stationary components
est_st <- Est_FDF(argvals=tt,
                  data=MatYield2,
                  stationary=TRUE,
                  h=19,k=2,p=5,
                  nbasis=15,
                  basis='Bspline',
                  kern_type = "BT",
                  lambda = NULL,
                  plot = FALSE)
# plots
plot(est_nst$hat.F, ylim = c(-2.5,2), main='F 1', xlab='s', ylab='')
plot(est_st$hat.F[1], ylim=c(-2.5,2), main='F 2', xlab='s', ylab='' )
plot(-est_st$hat.F[2], ylim=c(-2.5,2), main='F 3', xlab='s', ylab='')

# plot of time series
plot(1:N,est_nst$hat.beta, type='l', main='Factor process 1',
     xlab='n', ylab='')
plot(1:N,est_st$hat.beta[,1], type='l', main='Factor process 2',
     xlab='n', ylab='')
plot(1:N,est_st$hat.beta[,2], type='l', main='Factor process 3',
     xlab='n', ylab='')


