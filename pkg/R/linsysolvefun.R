##*************************************************************************
## linsysolvefun: Solve H*x = b
##
## x = linsysolvefun(L,b)
## where L contains the triangular factors of H. 
## 
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##*************************************************************************
 
linsysolvefun=function(L,b){
 
  if(is.null(dim(b))){
    x <- rep(0,length(b)) 
    if(L$matfct_options=="chol")
      x[L$perm] <- mextriang(L$L, mextriang(L$L,b[L$perm],2), 1) 
    else if(L$matfct_options=="spchol"){
      x[L$perm] <- mextriangsp(L$Lt,mextriangsp(L$L,b[L$perm],2), 1) 
    }else if(L$matfct_options=="lu")
      x <- solve(L$u)%*%(solve(L$l)%*%b[L$p@perm])
    else if(L$matfct_options=="splu"){     
      x[L$q@perm] <- (solve(L$u)%*%solve(L$l)%*%b[L$p@perm])[,1]
    }
  }else{
    x <- zeros(dim(b)[1],dim(b)[2]) 
    for(k in 1:(dim(b)[2])){
      if(L$matfct_options=="chol")
        x[L$perm,k] <- mextriang(L$L, mextriang(L$L,b[L$perm,k],2), 1) 
      else if(L$matfct_options=="spchol"){
        x[L$perm,k] <- mextriangsp(L$Lt,mextriangsp(L$L,b[L$perm,k],2), 1) 
      }else if(L$matfct_options=="lu")
        x[,k] <- solve(L$u)%*%(solve(L$l)%*%b[L$p@perm,k])
      else if(L$matfct_options=="splu"){     
        x[L$q@perm,k] <- (solve(L$u)%*%solve(L$l)%*%b[L$p@perm,k])[,1]
      }
    }
  }
  return(x)
}
