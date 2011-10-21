##*******************************************************************
## schurmat_qblk: compute schur matrix corresponding to SOCP blocks.
##
## NT  direction: output = schur + Ae*Ae' - Ad*Ad'
##
## where schur = A*D*A', and Ad is the modification to ADA' 
## so that the latter is positive definite. 
##
## [schur,UU,EE] = schurmat_qblk(blk,A,schur,UU,EE,p);
## 
## UU: stores the dense columns of Ax, Ae, Ad, and possibly 
##     those of A*D^{1/2}. It has the form UU = [Ax Ae Ad]. 
## EE: stores the assocaited (2,2) block matrix when the
##     output matrix is expressed as an augmented matrix.
##     It has the form EE = [0 -lam 0; -lam 0 0; 0 0 I].
## 
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##*******************************************************************

schurmat_qblk =   function(blk,A,par,schur,UU,EE,p){
	
	if(is.null(EE))
		count <- 0
	else 
		count <- max(max(EE[,2]),max(EE[,1])) 
	
	n <- sum(blk$size[[p]])
	numblk <- length(blk$size[[p]]) 
	blks <- list()
	blks$type <- blk$type[p]
	blks$size <- blk$size[[p]]
	
	ddsch <- par$dd[[p]]       
	Ae <- qprod(blks,A[[p]],par$ee[[p]]) 
	idxden <- checkdense(Ae)
	if(!is.null(idxden)){
		spcolidx <- setdiff(c(1:numblk),idxden) 
		s <- 1 + c(0,cumsum(blk$size[[p]]))
		idx <- s[idxden]
		tmp <- zeros(n,1) 
		tmp[idx] <- sqrt(2*abs(ddsch[idx])) 
		Ad <- qprod(blks,A[[p]],tmp)              
		ddsch[idx] <- abs(ddsch[idx])
		len <- length(idxden)
		w2 <- par$gamz[[p]]/par$gamx[[p]] 
		lam <- w2[idxden]
		if(len>1){
			spdiags <- as(diag(sqrt(lam)),"sparseMatrix")
			if(is.null(UU))
				UU <- cBind(Ae[,idxden]%*%spdiags,Ad[,idxden]) 
			else
				UU <- cbind(UU,Ae[,idxden]%*%spdiags,Ad[,idxden]) 
		}else{
			UU <- cbind(UU,Ae[,idxden]*sqrt(lam),Ad[,idxden])
		}
		tmp <- count+c(1:len) 
		EE <- rbind(EE, rbind(cbind(tmp, tmp, -lam),
						cbind(len+tmp, len+tmp, ones(len,1)))) 
		count <- count + 2*len
		Ae <- Ae[,spcolidx]      
		schur <- schur + Ae%*%tt(Ae)
	}else{
		schur <- schur + Ae%*%tt(Ae)
	}
	idxdenAq <- checkdense(A[[p]]) 
	if(!is.null(idxdenAq)){
		idxden <- idxdenAq 
		len <- length(idxden)
		if(len>1){
			spdiags <- as(diag(sqrt(abs(ddsch[idxden]))),"sparseMatrix") 
			Ad <- (A[[p]])[,idxden]%*%spdiags
		}else{
			Ad <- (A[[p]])[,idxden]*sqrt(abs(ddsch[idxden]))
		}
		if(is.null(UU))
			UU <- Ad
		else
			UU <- cBind(UU,Ad)
		tmp <- count+c(1:len) 
		EE <- rbind(EE, cbind(tmp,tmp,-sign(ddsch[idxden]))) 
		count <- count + len
		ddsch[idxden] <- zeros(len,1)
	}
	spdiags <- as(diag(ddsch),"sparseMatrix")
	schurtmp <- A[[p]]%*%spdiags%*%tt(A[[p]]) 
	schur <- schur + schurtmp
	
	return(list(schur=schur,UU=UU,EE=EE))
}
