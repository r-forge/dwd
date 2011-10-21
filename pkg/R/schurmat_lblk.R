##*******************************************************************
## schurmat_lblk: 
##
## [schur,UU,EE] <- schurmat_lblk(blk,A,par,schur,UU,EE,p);
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##*******************************************************************

schurmat_lblk = function(blk,A,par,schur,UU,EE,p){
  n <- sum(blk$size[[p]])  
  idxdenAl <- checkdense(A[[p]]) 
  ddsch <- par$dd[[p]] 
  if(!is.null(idxdenAl)){
    idxden <- idxdenAl
    len <- length(idxden)
    if(len>1){
      spdiags <- as(diag(sqrt(ddsch[idxden])),"sparseMatrix")
      Ad <- (A[[p]])[,idxden]%*%spdiags
    }else{
      Ad <- (A[[p]])[,idxden]*sqrt(ddsch[idxden])
    }
    UU <- cBind(UU, Ad)
    if(is.null(EE))
      count <- 0 
    else
      count <- max(max(EE[,1]),max(EE[,2])) 
    tmp <- count + c(1:len)
    EE <- rbind(EE,cbind(tmp,tmp,-ones(len,1))) 
    ddsch[idxden] <- zeros(len,1) 
  }
  spdiags <- as(diag(as.vector(ddsch[,1])),"sparseMatrix")
  schurtmp <- A[[p]]%*%spdiags%*%tt(A[[p]]) 
  schur <- schur + schurtmp
  return(list(schur=schur,UU=UU,EE=EE))
}
