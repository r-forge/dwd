##**********************************************************************
## NTscaling: Compute NT scaling matrix
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##**********************************************************************

NTscaling = function(blk,X,Z){

  numblk <- length(blk$type)
  gamx <- list()
  gamz <- list()
  dd <- list()
  ee <- list()
  ff <- list()

  for(p in 1:length(blk$type)){
    numblk <- length(blk$size[[p]])  
    n = sum(blk$size[[p]])  
    if(blk$type[p]=="l")
      dd[[p]] <- X[[p]]/Z[[p]]       
    else if(blk$type[p]=="q"){  
      blks <- list()
      blks$size <- blk$size[[p]]
      gamx[[p]] <- sqrt(qops(blks,X[[p]],X[[p]],2)) 
      gamz[[p]] <- sqrt(qops(blks,Z[[p]],Z[[p]],2)) 
      w2 <- gamz[[p]]/gamx[[p]]
      w <- sqrt(w2) 
      dd[[p]] <- qops(blks,1/w2,ones(n,1),4)
      tt <- qops(blks,1/w,Z[[p]],3) - qops(blks,w,X[[p]],4)
      mtt <- qops(blks,tt,tt,2)
      mtt[which(mtt<1e-16)] <- 1e-16
      gamtt <- sqrt(mtt) 
      ff[[p]] <- qops(blks,1/gamtt,tt,3) 
      ee[[p]] <- qops(blks,sqrt(2)/w,ff[[p]],4) 
    }
  }
  return(list(gamx=gamx,gamz=gamz,dd=dd,ee=ee,ff=ff))
}
