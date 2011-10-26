##**********************************************************************
## NTpred: Compute (dX,dy,dZ) for NT direction. 
##                       
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##**********************************************************************

NTpred = function(blk,A,rp,Rd,sigmu,X,Z,global_var){
	
	# compute NT scaling matrix
	
	par <- list()
	ntscalingf <- NTscaling(blk,X,Z)
	
	par$gamx <- ntscalingf$gamx
	par$gamz <- ntscalingf$gamz
	par$dd <- ntscalingf$dd
	par$ee <- ntscalingf$ee
	par$ff <- ntscalingf$ff
	
	# compute schur matrix
	m <- length(rp)
	schur <- as(matrix(0,m,m),"sparseMatrix") 
	UU <- NULL
	EE <- NULL
	Afree <- NULL 
	dX <- list()
	dy <- NULL
	dZ <- list()
	
	for(p in 1:length(blk$type)){
		if(blk$type[p]=="q"){       
			schurf <- schurmat_qblk(blk,A,par,schur,UU,EE,p)
			schur <- schurf$schur
			UU <- schurf$UU
			EE <- schurf$EE
		}else if(blk$type[p]=="l"){
			schurf <- schurmat_lblk(blk,A,par,schur,UU,EE,p)
			schur <- schurf$schur
			UU <- schurf$UU
			EE <- schurf$EE
		}else if(blk$type[p]=="u"){            
			if(ncol(A[[p]]==1))
				Afree = A[[p]]
			else
				Afree = cbind(Afree,A[[p]])
			
		}
	}
	
	# compute rhs
	
	ntrhsf <- NTrhsfun(blk,A,par,X,Z,rp,Rd,sigmu)
	rhs <- ntrhsf$rhs
	EinvRc <- ntrhsf$EinvRc
	hRd <- ntrhsf$hRd
	
	# solve linear system
	
	linsolvef <- linsysolve(schur,UU,Afree,EE,rhs,global_var)
	xx <- linsolvef$xx
	coeff <- linsolvef$coeff
	L <- linsolvef$L
        global_var = linsolvef$global_var

	# compute (dX,dZ)
	
	ntdirf <- NTdirfun(blk,A,par,Rd,EinvRc,xx,global_var)
	return(list(par=par,dX=ntdirf$dX,dy=ntdirf$dy,dZ=ntdirf$dZ,
                    coeff=coeff,L=L,hRd=hRd,global_var=ntdirf$global_var))
}
