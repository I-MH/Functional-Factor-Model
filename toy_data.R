
#This code simulates toy data from model (3) on the main paper

source('DataSimulation.R')

N <- 200  # number of curves to be simulated

#############################################################
# Stationary case:
set.seed(1579)
fn1 <- arima.sim(list(order=c(1,0,0),ar=0.6),n=N)
fn2 <- arima.sim(list(order=c(1,0,0),ar=-0.6 ), n=N)
ff <- cbind(fn1,fn2)
toy.data1 <- FDF.sim(m=50, rangeval=c(0,1), typew = 'MB',fn=ff)

par(mfrow=c(3,1))
matplot(toy.data1$tt, toy.data1$loadfMat, type = 'l', main='Factor loadings used to simulate',
        xlab = 's', ylab = '')
matplot(ff, type = 'l', main='Time series used to simulate',
        xlab = 'n')
matplot(toy.data1$tt, toy.data1$Data, type = 'l', lty=1, 
        main='Functional data simulated (nonstationary)',
        xlab = 's', ylab = '')
par(mfrow=c(1,1))

save(toy.data1, file="ToyData1.RData")

#############################################################
# Nonstationary case:
set.seed(1779)
fn1 <- arima.sim(list(order=c(1,1,0),ar=0.7 ), n=N-1) # I(1) process
fn2 <- arima.sim(list(order=c(1,0,0),ar=0.7 ), n=N)
ff <- cbind(fn1,fn2)

toy.data2 <- FDF.sim(m=50, rangeval=c(0,1), typew = 'MB',fn=ff)

par(mfrow=c(3,1))
matplot(toy.data2$tt, toy.data2$loadfMat, type = 'l', main='Factor loadings used to simulate',
        xlab = 's', ylab = '')
matplot(ff, type = 'l', main='Time series used to simulate',
        xlab = 'n')
matplot(toy.data2$tt, toy.data2$Data, type = 'l', lty=1, 
        main='Functional data simulated (nonstationary)',
        xlab = 's', ylab = '')
par(mfrow=c(1,1))
save(toy.data2, file="ToyData2.RData")


