##*********************************************************
## Atyfun: compute At*y
##
##  Q = Atyfun(blk,A,y);
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##**********************************************************

Atyfun=function(blk,A,y){
  Q <- list()
  for(p in 1:length(blk$type)){
    Q[[p]] <- A[[p]]%*%y 
  }
  return(Q)
}
