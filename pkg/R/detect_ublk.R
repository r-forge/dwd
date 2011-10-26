##*******************************************************************
## detect_ublk: search for implied free variables in linear
##              block. 
## [blk2,At2,C2,ublkinfo] = detect_ublk(blk,At,C); 
## 
## Important: blk is assumed to have only 1 linear block.
##
## i1,i2: indices corresponding to splitting of unrestricted varaibles
## i3   : remaining indices in the linear block
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##*******************************************************************

detect_ublk = function(blk,At,C){

  set.seed(0)   
  blk2 <- blk
  At2 <- At
  C2 <- C 
  numblk <- length(blk$type)
  ublkinfo <- list()
  for(p in 1:numblk){
    m <- dim(At[[p]])[2]        
    if(blk$type[p]=="l"){
      r <- rnorm(m)
      stime <- proc.time()
      Ap <- t(At[[p]])
      Cp <- C[[p]]
      ApTr <- t(r%*%Ap)
      sApTr <- sort(abs(ApTr))
      perm <- order(abs(ApTr))
      idx0 <- which(abs(diff(sApTr)) < 1e-14);
      i1 <- NULL
      i2 <- NULL
      if(!is.null(idx0)){
        n <- blk$size[[p]] 
        i1 <- perm[idx0]
        i2 <- perm[idx0+1]
        Api1 <- Ap[,i1]
        Api2 <- Ap[,i2]
        Cpi1 <- t(Cp[i1])
        Cpi2 <- t(Cp[i2])
        idxzr <- which(abs(Cpi1+Cpi2)<1e-14&apply(abs(Api1+Api2),2,sum)<1e-14)
        if(!is.null(idxzr)){
          i1 <- i1[idxzr]
          i2 <- i2[idxzr]
          blk2$type[p] <- "u"; 
          blk2$size[[p]] <- length(i1) 
          At2[[p]] <- t(Ap[,i1])
          C2[[p]]  <- Cp[i1]
          print(c("linear variables from unrestricted variable",
                  2*length(i1))) 
          i3 <- setdiff(c(1:n),union(i1,i2))
          if(!is.null(i3)){
            blk2$type[numblk+1] <- "l" 
            blk2$size[numblk+1] <- length(i3)
            At2[[numblk+1]] <- t(Ap[,i3])
            C2[[numblk+1]]  <- Cp[i3] 
          }
          ublkinfo$i1[[p]] <- i1
          ublkinfo$i2[[p]] <- i2
          ublkinfo$i3[[p]] <- i3
        }
      }
    }
  }
  return(list(blk2=blk2,At2=At2,C2=C2,ublkinfo=ublkinfo))
}
