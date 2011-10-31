##********************************************************************
## infeaspt: generate an initial point for sdp.m
##
##  [X0,y0,Z0] = infeaspt(blk,At,C,b,options,scalefac);
##
##  options = 1  if want X0,Z0 to be scaled identity matrices
##          = 2  if want X0,Z0 to be scalefac*(identity matrices).
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##********************************************************************

infeaspt=function(blk,At,C,b,options,scalefac,spdensity){

  options <- 1
  if (options == 1)
    scalefac <- NULL
  
  X0 <- list()
  Z0 <- list()
  m <- length(b) 
  
  myval <- validate(blk,At,C,b,spdensity=spdensity)
  At <- myval$At
  C <- myval$C
  spdensity = myval$spdensity
  
  for(p in 1:length(blk$type)){
    blktmp <- blk$size[[p]]
    n <- length(C[[p]])
    y0 <- zeros(m,1)
    b2 <- 1 + abs(t(b))
    if (options == 1){
      if(blk$type[p]=="q"){
        s <- 1+c(0,cumsum(blktmp))
        len <- length(blktmp)
        normC <- 1+normsvd(C[[p]])
        normA <- 1+sqrt(apply(At[[p]]*At[[p]],2,sum))
        idenqX <- zeros(sum(blktmp),1)
        idenqX[s[1:len]] <- sqrt(t(blktmp))
        X0[[p]] <- max(c(1,max(b2/normA)))*idenqX
        idenqZ <- zeros(sum(blktmp),1)
        normax <- max(c(normA,normC)) 
        idenqZ[s[1:len]] <- t(apply(rbind(sqrt(blktmp),normax*ones(1,len)),
                                    2,max))
        Z0[[p]] <- idenqZ
      }else if(blk$type[p]=="l"){
        normC <- 1+normsvd(C[[p]])
        normA <- 1+sqrt(apply(At[[p]]*At[[p]],2,sum))
        X0[[p]] <- max(1,max(b2/normA))*ones(n,1);
        Z0[[p]] <- max(1,max(c(normA,normC))/sqrt(n))*ones(n,1);
      }else if(blk$type[p]=="u"){
        X0[[p]] <- as(zeros(n,1),"sparseMatrix")
        Z0[[p]] <- as(zeros(n,1),"sparseMatrix")
      }else{
        print(" blk: some fields not specified correctly") 
      }
    }else if (options == 2){
      if(blk$type[p]=="q"){
        s <- 1+c(0,cumsum(blktmp))
        len <- length(blktmp)
        idenq <- zeros(sum(blktmp),1)
        idenq[s[1:len]] <- ones(len,1)
        X0[[p]] <- scalefac*idenq
        Z0[[p]] <- scalefac*idenq
      }else if(blk$type[p]=="l"){
        X0[[p]] <- scalefac*ones(n,1)
        Z0[[p]] <- scalefac*ones(n,1)
      }else if(blk$type[p]=="u"){
        X0[[p]] <- as(zeros(n,1),"sparseMatrix")
        Z0[[p]] <- as(zeros(n,1),"sparseMatrix")
      }else{
        print(" blk: some fields not specified correctly")
      }
    }
  }
  return(list(X0=X0,y0=y0,Z0=Z0))
}
