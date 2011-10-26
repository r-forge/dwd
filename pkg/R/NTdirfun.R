##*******************************************************************
## NTdirfun: compute (dX,dZ), given dy, for the NT direction.
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##*******************************************************************

NTdirfun=function(blk,A,par,Rd,EinvRc,xx,global_var){
  
  dX <- list()
  dZ <- list()
  dy <- NULL
  if(any(is.nan(xx)) | any(is.infinite(xx))){
    solve_ok <- 0
    global_var$solve_ok = solve_ok
    print("linsysolve: solution contains NaN or inf.")
    return
  }
  m <- dim(A[[1]])[1] 
  dy <- xx[1:m]
  count <- m
  
  for(p in 1:length(blk$type)){
    if(blk$type[p]=="l"){
      dZ[[p]] <- Rd[[p]] - mexMatvec(A[[p]],dy,1)     
      tmp <- par$dd[[p]]*dZ[[p]]
      dX[[p]] <- EinvRc[[p]] - tmp
    }else if(blk$type[p]=="q"){
      blks <- list()
      blks$type <- blk$type[p]  
      blks$size <- blk$size[[p]]  
      dZ[[p]] <- Rd[[p]] - mexMatvec(A[[p]],dy,1)  
      tmp <- par$dd[[p]]*dZ[[p]] +
        qops(blks,qops(blks,dZ[[p]],par$ee[[p]],1),par$ee[[p]],3) 
      dX[[p]] <- EinvRc[[p]] - tmp       
    }else if(blk$type[p]=="u"){ 
      n = sum(blk$size[[p]]) 
      dZ[[p]] <- zeros(n,1) 
      dX[[p]] <- matrix(xx[count+c(1:n)],nrow=n) 
      count <- count + n  
    }
  }
  return(list(dX=dX,dy=dy,dZ=dZ,global_var=global_var))
}
