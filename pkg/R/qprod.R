##***************************************************
## qprod: 
##       
## Input: A = [A1 A2 ... An]
##        x = [x1; x2; ...; xn]
## Output: [A1*x1 A2*x2 ... An*xn]
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##***************************************************

qprod = function(pblk,A,x){
  
  if(is(x,"sparseMatrix"))
    x <- as(x,"denseMatrix")
  n <- length(x)
  ii <- c(1:n)
  jj <- mexexpand(pblk$size,c(1:length(pblk$size))) 
  X <- spMatrix(n,max(jj),i=ii, j=jj, x=x)
  Ax <- A%*%X

  return(Ax)
}
