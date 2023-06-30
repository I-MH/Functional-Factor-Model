# temporal file
#----------------------------------------
require(ggplot2)
require(latex2exp)
source('FDFv2.R')
source('hatK.R')
source('SmoothData.R')
source('data.fd.R')
source('lrc.R')

# load data
load("ToyData1.RData")

data=toy.data1$Data; h=19;k=2;p=5;nbasis=25;basis='Fourier'
argvals=NULL;kmax=6; kern_type = "BT"
rangeval =c(0,1);nbasis=25;lambda = NULL

est1 <- FDF2(data=toy.data1$Data,
            h=h,k=2,p=5,nbasis=25,
            basis='Bspline',
            stationary = TRUE,
            type_lrc = 'lrcFc')

plot(est1$hat.F, lty=1)

lines(toy.data1$tt, -toy.data1$loadfMat[,1], col=1, lwd=2, lty=2)
lines(toy.data1$tt, toy.data1$loadfMat[,2], col=2, lwd=2, lty=2)

plot(est1$Varmx$harmonics)





