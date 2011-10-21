##******************************************************************
## blkcholfun: compute Cholesky factorization of X. 
##          
##  [Xchol,indef] = blkcholfun(blk,X); 
##  
##  X = Xchol'*Xchol;
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##******************************************************************
 
blkcholfun = function(blk,X){
  if(!(class(X)=="list")){
    indef <- 0 
    n <- max(dim(X))
    if(blk$type=="q"){ 
      gamx <- mexqops(blk,X,X,2) 
      if(any(gamx <= 0)) indef <- 1 
      Xchol <- NULL   
    }else if(blk$type=="l"){ 
      if(any(X <= 0)) indef <- 1
      Xchol <- NULL  
    }else if(blk$type=="u"){
      Xchol <- NULL
    }
  }else{ 
    Xchol <- list()
    indef <- NULL
    for(p in 1:length(blk$type)){
      blklist <- list(type=blk$type[p],size=blk$size[[p]])
      cholfun <- blkcholfun(blklist,X[[p]])
      Xchol[[p]] <- cholfun$Xchol
      indef[p] <- cholfun$indef
    }
    indef <- max(indef)
  }
  return(list(Xchol=Xchol,indef=indef))
}
