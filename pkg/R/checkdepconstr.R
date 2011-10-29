##*****************************************************************************
## checkconst: compute AAt to determine if the 
##             constraint matrices Ak are linearly independent. 
##              
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##*****************************************************************************

checkdepconstr =  function(blk,At,b,y,rmdepconstr){
   
  ## compute AAt
  m <- length(b)
  AAt <- as(zeros(m,m),"sparseMatrix")
  numdencol <- 0
  UU <- NULL
  for(p in 1:length(blk$type)){
    decolidx <- checkdense(tt((At[[p]]))) 
    if(!is.null(decolidx)){ 
      n2 <- dim(At[[p]])[1] 
      dd <- ones(n2,1) 
      len<- length(decolidx); 
      dd[decolidx] <- zeros(len,1)
      AAt <- AAt + tt((At[[p]]))%*%Diagonal(n2,dd)%*%At[[p]]
      tmp <- tt((At[[p]])[decolidx,])
      UU <- cbind(UU,tmp) 
      numdencol <- numdencol + len 
    }else{
      AAt <- AAt + tt((At[[p]]))%*%At[[p]]  
    }
  }
  
  if (numdencol > 0)
    print(c("number of dense column in AAt = ",numdencol)) 

  numdencol <- dim(UU)[2]

  if(!is(AAt,"sparseMatrix"))
    AAt <- as(AAt,"sparseMatrix")
  pertdiag <- 1e-16*normsvd(AAt)*ones(m,1)
  AAt <- AAt + Diagonal(m,pertdiag)
  mychol <- Cholesky(AAt)
  if (class(mychol)=="try-error"){
    print(" AAt is not positive definite ")
    indef <- -1
    return
  }else{
    indef <- 0
    r <- expand(mychol)
    perm <- r$P@perm
    R <- tt(r$L[perm,])
    dd <- Diag(r$L)       
  }

  ## find independent rows of A

  idxB <- c(1:m)   
  feasible <- 1
  ndepconstr <- 0
  idxN <- which(dd[,1] < 1e-8)
  if(length(idxN)!=0){
    print(c("number of nearly dependent constraints = ",length(idxN)))
    ndepconstr <- 1
    if (numdencol==0){
      if (rmdepconstr){
        idxB <- setdiff(c(1:m),idxN)     
        print("checkdepconstr: removing dependent constraints...")
        W <- findcoeff(blk,At,idxB,idxN)
        tmp <- tt(W)%*%b[idxB] - b[idxN];
        nnorm <- normsvd(tmp)/max(1,normsvd(b)) 
        tol <- 1e-8
        if (nnorm > tol){
          feasible <- 0 
          print("checkdepconstr: inconsistent constraints exist,")
          print("problem is infeasible.");
        }else{
          feasible <- 1;
          for(p in 1:length(blk)){ 
            At[[p]] <- (At[[p]])[,idxB]
          }
          b <- b[idxB]
          y <- y[idxB]
        }
      }else{
        print("To remove these constraints,")
        print(" re-run sqlp.m with OPTIONS.rmdepconstr = 1")
      }
    }else{
      print("warning: the sparse part of AAt may be nearly");
      print("singular.");
    }
  }
  return(list(At=At,b=b,y=y,idxB=idxB,
              ndepconstr=ndepconstr,feasible=feasible))
}

##
##*****************************************************************************
## findcoeff: 
##
## W = findcoeff(blk,At,idXB,idXN);
## 
## idXB = indices of independent columns of At. 
## idxN = indices of   dependent columns of At.
## 
## AB = At(:,idxB); AN = At(:,idxN) = AB*W
##
## 
## SDPT3: version 3.0
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last modified: 2 Feb 01
##*****************************************************************************

findcoeff=function(blk,At,idxB,idxN){

  tol <- 1e-8; 
  AB <- NULL
  AN <- NULL
  for(p in 1:length(blk$size)){
    AB <- rbind(AB,(At[[p]])[,idxB])
    AN <- rbind(AN,(At[[p]])[,idxN])
  }
  m <- dim(AB)[1]
  n <- dim(AB)[2]

  ##-----------------------------------------
  ## find W so that AN = AB*W
  ##-----------------------------------------

  elu <- expand(lu(as(AB,"sparseMatrix")))
  L <- elu$L
  U <- elu$U
  P <- elu$P
  Q <- elu$Q
  rhs  <- P%*%AN
  Lhat <- L[1:n,] 
  W <- Q%*%solve(U)%*%solve(Lhat)%*%rhs[1:n,]
  nnorm <- normsvd(AN-AB%*%W)/max(1,normsvd(AN))
  if (nnorm > tol) 
    print("\n findcoeff: basis rows may be identified incorrectly.")
  return(W)

}
tt <- function(x){
  if(is(x,"sparseMatrix"))
    return(as(t(as.matrix(x)),"sparseMatrix"))
  else if(is(x,"gdeMatrix"))
    return(t(as.matrix(x)))
  else
    return(t(x))
}
Diag <- function(x){
  if(is(x,"sparseMatrix"))
    return(as(diag(as.matrix(x)),"sparseMatrix"))
  else
    return(diag(x))
}

