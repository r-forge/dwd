##*********************************************************
## AXfun: compute A*X
##
##   AX = AXfun(blk,A,X);
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##**********************************************************

AXfun = function(blk,A,X){
  m <- dim(A[[1]])[1] 
  AX <- zeros(m,1)
  for(p in 1:length(blk$type)){
    AX <- AX + A[[p]]%*%X[[p]]
  }
  return(AX)
}
