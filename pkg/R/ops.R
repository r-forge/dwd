##******************************************************************
## ops:  
##
##   Z = ops(X,operand,Y,alpha); 
##
##  INPUT:        X = a matrix or a scalar
##                    or a CELL ARRAY consisting only of matrices 
##          operand = sym, transpose, triu, tril,
##                    real, imag, sqrt, abs, max, min, nnz,
##                    spdiags, ones, norm, sum, row-norm, 
##                    rank1, rank1inv, inv
##                    +,  -, *, .*,  ./, .^ 
##     Y = empty or a matrix or a scalar 
##               or a CELL ARRAY consisting only of matrices
##    alpha  = empty or a scalar 
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##******************************************************************

ops=function(X,operand,Y,alpha){

  if (nargs()==2){
    if(operand=="norm"){
      Z <- 0
      for(p in 1:length(X))
        Z = Z + sum(sum(X[[p]]*X[[p]]))
      Z <- sqrt(Z)
    }else if(operand=="getM"){
      Z <- 0
      for(p in 1:length(X))
        Z <- Z + dim(X[[p]])[1]
    }
  }else if (nargs()==3){
    if ((operand=="+")|(operand=="-")|(operand==".*")){ 
      Z <- list() 
      if(operand=="+"){
        for(p in 1:length(X))
          Z[[p]] <- X[[p]] + Y[[p]]
      }else if(operand=="-"){
        for(p in 1:length(X))
          Z[[p]] <- X[[p]] - Y[[p]]
      }else if(operand==".*"){
        for(p in 1:length(X))
          Z[[p]] <- X[[p]]* Y[[p]]
      }else{
        print(c(operand,"is not available"))
      }
    }
  }else{
    Z = list() 
    if(operand=="+"){ 
      for(p in 1:length(X))
        Z[[p]] <- X[[p]] + alpha*Y[[p]]
    }else if(operand=="-"){ 
      for(p in 1:length(X))
        Z[[p]] <- X[[p]] - alpha*Y[[p]]
    }else{
      print(c(operand,"is not available"))
    }
  }
  return(Z)
}
