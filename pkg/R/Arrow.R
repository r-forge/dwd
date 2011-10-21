##********************************************************
## Arrow: 
##  
## Fx = Arrow(pblk,f,x,options); 
##
## if options == 0; 
##    Fx = Arr(F)*x
## if options == 1; 
##    Fx = Arr(F)^{-1}*x 
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##********************************************************

Arrow=function(blk,f,x,options){
  if(nargs()==3)
    options <- 0
  s <- 1 + c(0,cumsum(blk$size)) 
  idx1 <- s[1:length(blk$size)]
  if(options==0){
    inprod <- mexqops(blk,f,x,1)  
    Fx <- mexqops(blk,f[idx1],x,3) + mexqops(blk,x[idx1],f,3) 
    Fx[idx1] <- inprod 
  }else{
    gamf2 <- mexqops(blk,f,f,2)
    gamprod <- mexqops(blk,f,x,2)
    alpha <- gamprod/gamf2 
    Fx <- mexqops(blk,1/f[idx1],x,3) - mexqops(blk,alpha/f[idx1],f,3) 
    Fx[idx1] <- alpha
  }
  return(Fx)
}
