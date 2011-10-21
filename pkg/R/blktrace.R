##**********************************************************************
## blktrace: compute <X1,Z1> + ... <Xp,Zp>
##              
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##**********************************************************************

blktrace=function(blk,X,Z){
  trXZ <- 0 
  for(p in 1:length(blk$type)){
    trXZ <- trXZ + tt(as.matrix(X[[p]]))%*%Z[[p]] 
  }
  return(trXZ)
}
