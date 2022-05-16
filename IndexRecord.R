
#-------------------------------------------------------------------------------
# functions from R packages 'fda' and 'fdaoutlier'. These are references on 
# the main paper

combinat=function(n,p){
  if (n<p){combinat=0}
  else {combinat=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}
}
#BD2
fBD2=function(data){
  p=dim(data)[1]
  n=dim(data)[2]
  rmat=apply(data,1,rank)
  down=apply(rmat,1,min)-1
  up=n-apply(rmat,1,max)
  (up*down+n-1)/combinat(n,2)
}
#MBD
fMBD=function(data){
  p=dim(data)[1]
  n=dim(data)[2]
  rmat=apply(data,1,rank)
  down=rmat-1
  up=n-rmat
  (rowSums(up*down)/p+n-1)/combinat(n,2)
}

extremal_depth <- function(dts){
  if(is.data.frame(dts)){
    dts <- as.matrix(dts)
  }
  
  if(!is.array(dts) || !is.numeric(dts))
    stop("Argument \"dts\" must be a numeric matrix or dataframe.")
  
  if (any(!is.finite(dts))){
    stop("Missing or infinite values are not allowed in argument \"data\"")
  }
  if(nrow(dts) < 3) stop("The number of curves must be greater than 2")
  
  ddim <- dim(dts)
  n <- ddim[1]
  p <- ddim[2]
  pwdepth <- pwise_depth(dt = dts, n = n) # matrix of n by p
  pmfs <- apply(pwdepth, 1, function(x){
    pmf <- table(x)/p
    return(c(as.numeric(names(pmf[1])), # depth level and mass
             pmf[1]))
  })
  
  depth_levels <- pmfs[1,]
  masses <- pmfs[2, ]
  # order functions according to depth_levels and mass
  ordered_functions <- sapply(sort(unique(depth_levels),
                                   method = "quick"),
                              function(x){
                                fns_depth_level <- which(depth_levels == x)
                                if(length(fns_depth_level) > 1){
                                  fns_depth_level[order(masses[fns_depth_level], decreasing = T)]
                                }else{
                                  fns_depth_level
                                }
                              })
  ordered_functions <- unlist(ordered_functions)
  depth_values <- ((1:n)/n)[order(ordered_functions)]
}

pwise_depth <- function(dt, n) {
  pdepth <- apply(dt, 2, function(i){
    (1 - abs(2*rank(i) - n -1)/n) # for rank r, d = ((r - 1) + (n - r))/n
  })
  return(pdepth)
}

#-------------------------------------------------------------------------------


IndexRecord <- function(data, depth=c('MBD','ED'), plot=TRUE){
  # It computes the index time where record curves are observed
  # args:
  #   data: a matrix with dim m x n, where m is the total num of points observed for a curve
  #         and n is the total number of curves
  #   depth: depth to be use to order curves and estimate the records
  #
  # values: list 
  #   Est: a estatistic value 
  #   Records: Index/time where records are estimated
  #   UpperR: Index/time where upper record curves are estimated
  #   LowerR: Index/time where lower record curves  are estimated
  method <- match.arg(depth)
  N <- dim(data)[2]
  p <- dim(data)[1]
  Ind.Record <- rep(0,N)
  Ind.RecordU <- rep(0,N)
  Ind.RecordL <- rep(0,N)
  Ind.Record[1:2] <- 1 # the the first two are record curves 
  tlastR <- 2          # time corresponding to the last record
  indexminmax <- matrix(c(1,2),2,2)  # a matrix containing the two most extreme curves at each time
  Mat_means0 <-  rowMeans(data[,1:2])
  Ind.RecordL[1] <- ifelse( sum( data[,1]<Mat_means0 ) < floor(p/2), 0,1)
  Ind.RecordL[2] <- ifelse( sum( data[,2]<Mat_means0) < floor(p/2), 0,1)
  
  if(method=='MBD'){
    for (k in 3:N) {
      new.data=data[,1:k]
      depth=fMBD(new.data)
      index=order(depth)
      if( depth[k]<depth[tlastR] || depth[k]<depth[indexminmax[1,k-1]] || depth[k]<depth[indexminmax[2,k-1]]){
        Ind.Record[k] <- 1
        Ind.RecordL[k] <-  ifelse( sum( data[,k]<rowMeans(data[,1:k]) )< floor(p/2), 0,1)
        tlastR <- k
      }
      indexminmax <- cbind(indexminmax, index[1:2])
    }
  }else{
      for (k in 3:N) {
        new.data=data[,1:k]
        depth= extremal_depth(t(new.data))
        index=order(depth)
        if( depth[k]<depth[tlastR] || depth[k]<depth[indexminmax[1,k-1]] || depth[k]<depth[indexminmax[2,k-1]]){
          Ind.Record[k] <- 1
          Ind.RecordL[k] <-  ifelse( sum( data[,k]<rowMeans(data[,1:k]) )< floor(p/2), 0,1)
          tlastR <- k
        }
        indexminmax <- cbind(indexminmax, index[1:2])
      }
  }
  
  Ind.RecordU <- Ind.Record -Ind.RecordL 
  timeC <- 1/sqrt(1:N)
  Ind.S=timeC*cumsum(Ind.Record)
  
  if(plot){
    par(mfrow=c(1,2))
    plot(1:length(Ind.Record), cumsum(Ind.Record), type = 's', col=1,
         ylab = 'num of records', xlab = 'index of time', main=paste('Record times with', method))
    lines(1:length(Ind.RecordL), cumsum(Ind.RecordL), type = 's', col=4 )
    lines(1:length(Ind.RecordU), cumsum(Ind.RecordU), type = 's', col=2 )
    legend("topleft", legend = c('total records','upper records', "lower records" ), col = c(col=1,2,4),
           ncol = 1, cex = 1, lwd = 2, bty='n')
    matplot(data, type = 'l', col = grey(.7,.4), main='Functional Data', xlab = 's')
    matplot(data[,Ind.RecordU==1], type = 'l', col = 2, add = TRUE )
    matplot(data[,Ind.RecordL==1], type = 'l', col = 4, add = TRUE )
    par(mfrow=c(1,1))
  }
  
  return(list(Est=Ind.S, Records=Ind.Record, UpperR=Ind.RecordU, 
              LowerR=Ind.RecordL,IndexSupInf=indexminmax) )
}

