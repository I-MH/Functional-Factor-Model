
require(ggplot2)
require(latex2exp)
source('FDF.R')
source('hatK.R')
source('SmoothData.R')
source('data.fd.R')
source('lrc.R')

# stationary case --------------------------------------------------------------

# load data
load("ToyData1.RData")
matplot(toy.data1$tt, toy.data1$Data, type = 'l', main='Functional data',
        xlab = 's', ylab = ' ')

N=dim(toy.data1$Data)[2]
# selecting the optimal window parameter
# this can be skipped
#require(fChange)
#h_datafd <- SmoothData(data = toy.data1$Data, type_basis='Fourier', 
#                       nbasis=15,lambda = NULL)$fd
#opt_bandwidth(h_datafd, "PR", "BT", is_change = FALSE)$hat_h_opt
# this is 9
# end of this can be skipped

# estimate F_k and beta_{ki} using discrete data.
est1 <- FDF(data=toy.data1$Data,
                h=9,k=2,p=5,nbasis=25,
                basis='Fourier',
                stationary = TRUE)

data.f <- data.frame(s=toy.data1$tt, 
                     f1=toy.data1$loadfMat[,1], 
                     f2=toy.data1$loadfMat[,2])
data.f$hatf1 <- eval.fd(toy.data1$tt, est1$hat.F[1])
data.f$hatf2 <- eval.fd(toy.data1$tt, est1$hat.F[2])

p <- ggplot(data = data.f, aes(x=s))
p+geom_line(aes(y=f1, color='Factor 1', linetype = "Truth"),size=1) + 
  geom_line(aes(y=f2, color='Factor 2',linetype = "Truth"),size=1) +
  geom_line(aes(y=hatf1, color='Factor 1',linetype = "Estimated")) +
  geom_line(aes(y=hatf2, color='Factor 2',linetype = "Estimated"))+
  annotate('text', x=c(.2,.46), y=c(.9,.9), size=5,
           label=c(TeX('$F_2 (s)$'),TeX('$F_1 (s)$') ) )+
  labs(y=' ', x = "s")+
  ggtitle(TeX('Estimators of Factor Loadings $F_k (s)$') )+  
  guides(color = FALSE) +
  scale_linetype_discrete(name = element_blank() )+
  theme(legend.position='top',
        legend.direction = "horizontal", 
        legend.text = element_text(size = 15),
        plot.title = element_text(lineheight =.8,face = "bold",size = 19,hjust = 0.5),
        axis.text=element_text(size=13),
        axis.title=element_text(size=17), 
        axis.title.y = element_text(angle = 0,  vjust = .5))

data2.f <- data.frame(time=1:N,
                      f1=toy.data1$fn[,1],
                      f2=toy.data1$fn[,2],
                      hatf1=est1$hat.beta[,1], 
                      hatf2=est1$hat.beta[,2])

p <- ggplot(data = data2.f, aes(x=time))
p+geom_line(aes(y=f1, color='Truth', linetype = "Truth"), size=.8) +
  geom_line(aes(y=hatf1, color='Estimated',linetype = "Estimated"), size=.6)+
  labs(y=' ', x = "n")+
  ggtitle(TeX('Factor Process $\\beta_{n,1}$') )+  
  guides(color = FALSE) +
  scale_linetype_discrete(name = element_blank() )+ 
  theme(legend.position='top',legend.direction = "horizontal",
        legend.text = element_text(size = 15),
        plot.title = element_text(lineheight =.8,face = "bold",size = 19,hjust = 0.5 ),
        axis.text=element_text(size=13),
        axis.title=element_text(size=17), 
        axis.title.y = element_text(angle = 0,  vjust = .5))

p+geom_line(aes(y=f2, color='Truth', linetype = "Truth"), size=.8) +
  geom_line(aes(y=hatf2, color='Estimated',linetype = "Estimated"), size=.6)+
  labs(y=' ', x = "n")+
  ggtitle(TeX('Factor Process $\\beta_{n,2}$') )+  
  guides(color = FALSE) +
  scale_linetype_discrete(name = element_blank() )+ 
  theme(legend.position='top',legend.direction = "horizontal",
        legend.text = element_text(size = 15),
        plot.title = element_text(lineheight =.8,face = "bold",size = 19,hjust = 0.5 ),
        axis.text=element_text(size=13),
        axis.title=element_text(size=17), 
        axis.title.y = element_text(angle = 0,  vjust = .5))


# non stationary case ----------------------------------------------------------

# load data
load("ToyData2.RData")
N=dim(toy.data2$Data)[2]

# selecting the optimal window parameter
# this can be skipped -----------------------------------
#require(fChange)
#h1_datafd <- SmoothData(data= t(apply(toy.data2$Data,1, diff)), 
#                        type_basis='Fourier', 
#                        nbasis=25,
#                        lambda = NULL)$fd
#opt_bandwidth(h1_datafd,  kern_type='PR',kern_type_ini='BT',is_change = FALSE)$hat_h_opt 

# this is 12
# end of this can be skipped ----------------------------

# estimate the parameters for the non stationary components
est2_ns <- FDF(data=toy.data2$Data,
                   h=12,k=1,p=4,nbasis=25,
                   basis='Fourier',
                   stationary = FALSE)

data.f <- data.frame(s=toy.data2$tt,
                     f1=toy.data2$loadfMat[,1], 
                     f2=toy.data2$loadfMat[,2])
data.f$hatf1 <- eval.fd(toy.data2$tt, est2_ns$hat.F[1])

# subtracts the non-stationary part 
Data2 <- toy.data2$Data - est2_ns$Xhat

# selecting the optimal window parameter for stationary component -
# this can be skipped
h2_datafd <- SmoothData(data= Data2, 
                        type_basis='Fourier', 
                        nbasis=25,
                        lambda = NULL)$fd
opt_bandwidth(h2_datafd,  kern_type='PR',kern_type_ini='BT',is_change = FALSE)$hat_h_opt 
# this is 12
# end of this can be skipped
#-----------------------------------------------------------------

est2_s <-  FDF(data=Data2,
                   h=12,k=1,p=3,
                   nbasis=25,
                   basis='Fourier',
                   stationary = TRUE)

data.f$hatf2 <- eval.fd(toy.data2$tt, est2_s$hat.F[1])

p <- ggplot(data = data.f, aes(x=s))
p+geom_line(aes(y=f1, color='Factor 1', linetype = "Truth"), size=1) + 
  geom_line(aes(y=f2, color='Factor 2',linetype = "Truth"),size=1) +
  geom_line(aes(y=hatf1, color='Factor 1',linetype = "Estimated")) +
  geom_line(aes(y=hatf2, color='Factor 2',linetype = "Estimated"))+
  annotate('text', x=c(.2,.46), y=c(.9,.9), size=5,
           label=c(TeX('$F_2 (s)$'),TeX('$F_1 (s)$') ) )+
  labs(y=' ', x = "s")+
  ggtitle(TeX('Estimators of Factor Loadings $F_k (s)$') )+  
  guides(color = FALSE) +
  scale_linetype_discrete(name = element_blank() )+
  theme(legend.position='top',legend.direction = "horizontal",
        legend.text = element_text(size = 15),
        plot.title = element_text(lineheight =.8,face = "bold",size = 19,hjust = 0.5 ),
        axis.text=element_text(size=13),
        axis.title=element_text(size=17), 
        axis.title.y = element_text(angle = 0,  vjust = .5))

data2.f <- data.frame(time=1:N,f1=toy.data2$fn[,1],f2=toy.data2$fn[,2],
                      hatf1=est2_ns$hat.beta[,1], hatf2=est2_s$hat.beta[,1])
p <- ggplot(data = data2.f, aes(x=time))
p+geom_line(aes(y=f1, color='Truth', linetype = "Truth"), size=.8) +
  geom_line(aes(y=hatf1, color='Estimated',linetype = "Estimated"), size=.6)+
  labs(y=' ', x = "n")+
  ggtitle(TeX('Factor Process $\\beta_{n,1}$') )+  
  guides(color = FALSE) +
  scale_linetype_discrete(name = element_blank() )+ 
  theme(legend.position='top',legend.direction = "horizontal",
        legend.text = element_text(size = 15),
        plot.title = element_text(lineheight =.8,face = "bold",size = 19,hjust = 0.5 ),
        axis.text=element_text(size=13),
        axis.title=element_text(size=17), 
        axis.title.y = element_text(angle = 0,  vjust = .5))

p+geom_line(aes(y=f2, color='Truth', linetype = "Truth"), size=.8) +
  geom_line(aes(y=hatf2, color='Estimated',linetype = "Estimated"), size=.6)+
  labs(y=' ', x = "n")+
  ggtitle(TeX('Factor Process $\\beta_{n,2}$') )+  
  guides(color = FALSE) +
  scale_linetype_discrete(name = element_blank() )+ 
  theme(legend.position='top',legend.direction = "horizontal",
        legend.text = element_text(size = 15),
        plot.title = element_text(lineheight =.8,face = "bold",size = 19,hjust = 0.5 ),
        axis.text=element_text(size=13),
        axis.title=element_text(size=17), 
        axis.title.y = element_text(angle = 0,  vjust = .5))

