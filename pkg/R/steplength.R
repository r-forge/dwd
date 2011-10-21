##***************************************************************************
## steplength: compute xstep such that  X + xstep*dX >= 0.
##
## [xstep] = steplength(blk,X,dX);
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##***************************************************************************

steplength = function(blk,X,dX){
	eps <- 2.2204e-16
	xstep <- NULL
	for(p in 1:length(blk$type)){
		numblk <- length(blk$size[[p]]) 
		if(sum(blk$size[[p]]>1)){
			if((any(is.nan(dX[[p]][,1]))|any(is.infinite(dX[[p]][,1])))){
				xstep <- 0
				break
			}
		}else{
			if((is.nan(dX[[p]]))|(is.infinite(dX[[p]]))){
				xstep <- 0
				break
			}
		}
		if(blk$type[p]=="q"){
			blks <- list()
			blks$type <- blk$type[p]
			blks$size <- blk$size[[p]]
			aa <- qops(blks,dX[[p]],dX[[p]],2) 
			bb <- qops(blks,dX[[p]],X[[p]],2) 
			cc <- qops(blks,X[[p]],X[[p]],2)
			dd <- bb*bb - aa*cc 
			tmp <- apply(cbind(aa,bb),1,min)
			idx <- which(dd > 0 & tmp < 0)
			steptmp <- 1e12*ones(numblk,1) 
			if(!is.null(idx)){
				steptmp[idx] <- -(bb[idx]+sqrt(dd[idx]))/aa[idx]       
			}
			idx <- which(abs(aa) < eps & bb < 0) 
			if(!is.null(idx)){
				steptmp[idx] <- -cc[idx]/(2*bb[idx]) 
			}
			
			## also need first component to be non-negative
			
			ss <- 1 + c(0, cumsum(blk$size[[p]]))
			ss <- ss[1:length(blk$size[[p]])] 
			dX0 <- (dX[[p]])[ss]
			X0 <- (X[[p]])[ss]
			idx <- which(dX0 < 0 & X0 > 0)
			if(!is.null(idx)){
				steptmp[idx] <- apply(cbind(steptmp[idx],-X0[idx]/dX0[idx]),1,min)
			}
			xstep[p] <- min(steptmp)
		}else if(blk$type[p]=="l"){
			idx <- which(dX[[p]][,1] < 0) 
			if(!is.null(idx)&length(idx)>0){
				xstep[p] <- min(-X[[p]][idx]/dX[[p]][idx])
			}else{ 
				xstep[p] <- 1e12
			}
		}else if(blk$type[p]=="u"){
			xstep[p] <- 1e12
		}
	}
	xstep <- min(xstep) 
	return(xstep)
}
