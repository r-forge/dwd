##*******************************************************************
## NTrhsfun: compute the right-hand side vector of the 
##           Schur complement equation for the NT direction. 
## 
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##*******************************************************************

NTrhsfun=function(blk,A,par,X,Z,rp,Rd,sigmu,hRd,dX,dZ){
	
	m <- length(rp)   
	if (nargs() > 8) 
		corrector <- 1 
	else{ 
		corrector <- 0 
		hRd <- zeros(m,1) 
	}
	hEinvRc <- zeros(m,1) 
	rhsfree <- NULL
	EinvRc <- list() 
	for(p in 1:length(blk$type)){
		blks <- list()
		blks$type <- blk$type[p]
		blks$size <- blk$size[[p]]
		n <- sum(blk$size[[p]])
		numblk <- length(blk$size[[p]])  
		if (blk$type[p]=="l"){
			if (corrector)
				Rq <- dX[[p]]*dZ[[p]] 
			else{
				Rq <- as(zeros(n,1),"sparseMatrix") 
				tmp <- par$dd[[p]]*Rd[[p]]
				tmp2 <- mexMatvec(A[[p]],tmp,0) 
				hRd <- hRd + tmp2
			}
			EinvRc[[p]] <- sigmu/Z[[p]]-X[[p]] - Rq/Z[[p]]
			tmp2 <- mexMatvec(A[[p]],EinvRc[[p]],0)  
			hEinvRc <- hEinvRc + tmp2
		}else if(blk$type[p]=="q"){ 
			w <- sqrt(par$gamz[[p]]/par$gamx[[p]]) 
			if(corrector){
				hdx <- qops(blks,w,par$ff[[p]],5,dX[[p]]) 
				hdz <- qops(blks,w,par$ff[[p]],6,dZ[[p]]) 
#				nargin <<- 3
				hdxdz <- Arrow(blks,hdx,hdz)
				vv <- qops(blks,w,par$ff[[p]],5,X[[p]]) 
#				nargin <<- 4
				Vihdxdz <- Arrow(blks,vv,hdxdz,1) 
				Rq <- qops(blks,w,par$ff[[p]],6,Vihdxdz) 
			}else{
				Rq <- as(zeros(n,1),"sparseMatrix") 
				tmp <- par$dd[[p]]*Rd[[p]] +
						qops(blks,qops(blks,Rd[[p]],par$ee[[p]],1),
								par$ee[[p]],3)
				tmp2 <- mexMatvec(A[[p]],tmp,0)
				hRd <- hRd + tmp2
			}
			EinvRc[[p]] <- qops(blks,-sigmu/(par$gamz[[p]])^2,
					Z[[p]],4)-X[[p]]-Rq
			tmp2 <- mexMatvec(A[[p]],EinvRc[[p]],0)
			hEinvRc <- hEinvRc + tmp2
		}else if(blk$type[p]=="u"){
			if(ncol(Rd[[p]])==1)
				rhsfree = Rd[[p]]
			else
				rhsfree <- rbind(rhsfree,Rd[[p]]) 
		}
	}
	rhs <- rp + hRd - hEinvRc 
	if(!is.null(rhsfree)){
		if((ncol(rhsfree)==1)){
			rhs = rbind(as.matrix(rhs),rhsfree[1,1])
		}else
			rhs <- rbind(as.matrix(rhs), rhsfree)
	}else
		rhs <- as.matrix(rhs)
	return(list(rhs=rhs,EinvRc=EinvRc,hRd=hRd))
}
