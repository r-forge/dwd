##***************************************************************
## linsysolve: solve linear system to get dy, and direction
##             corresponding to unrestricted variables. 
##
## [xx,coeff,L,M] = linsysolve(schur,UU,Afree,EE,rhs); 
##
## child functions: symqmr.m, linsysolvefun.m
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##***************************************************************

linsysolve = function(schur,UU,Afree,EE,rhs,global_var){
	
        iter <- global_var$iter
        use_LU <- global_var$use_LU
        printlevel <- global_var$printlevel
        solve_ok <- global_var$solve_ok
        nnzmatold <- global_var$nnzmatold
        depconstr <- global_var$depconstr
        spdensity <- global_var$spdensity
        matfct_options <- global_var$matfct_options
        matfct_options_old <- global_var$matfct_options_old
	nnzmat <- 0
	if(iter==1) use_LU <- 0
	if(is.null(nnzmatold)) nnzmatold <- 0
	
	## diagonal perturbation, diagonally scale schur
	
	m <- max(dim(schur))
	diagschur <- abs(diag(as.matrix(schur)))
	dimafree <- 0
	if(!is.null(Afree))
		dimafree <- max(dim(Afree))
	if(dimafree>=0){
		if(depconstr)
			diagschur[which(diagschur<1)] <- 1
		else
			diagschur[which(diagschur<1e-4)] <- 1e-4      
		diagschur <- 1e-6*diagschur
		pertdiag <- as(diag(diagschur),"denseMatrix")
		#mexschurfun(schur,pertdiag,2) 
		schur <- schur + pertdiag
	}
	
	## assemble coefficient matrix
	
	len <- 0
	if(!is.null(Afree))
		len <- dim(Afree)[2]
	if(!is.null(EE))
		EE[,1:2] <- len + EE[,1:2] ## adjust for ublk
	if(len)
		EE = rbind(cbind(c(1:len),c(1:len),zeros(len,1)),EE) 
	
	coeff <- list()
	if(is.null(EE))
		coeff$mat22 <- NULL 
	else
		coeff$mat22 <- spMatrix(dim(EE)[1],dim(EE)[1],i=EE[,1],j=EE[,2],x=EE[,3])
	if(is.null(Afree))
		coeff$mat12 <- UU
	else{
		if(is.null(UU))
			coeff$mat12 <- Afree
		else
			coeff$mat12 <- cBind(Afree,UU)
	}
	coeff$mat11 <- schur ## important to use perturbed schur matrix
	ncolU <- 0
	if(!is.null(coeff$mat12)){
		ncolU <- dim(coeff$mat12)[2] 
		coeff$mat12 <- as(coeff$mat12,"dgTMatrix")
	}
	
	## pad rhs with zero vector
	## decide which solution methods to use
	
	rhs = c(rhs,zeros(m+ncolU-length(rhs),1)) 
	if (ncolU > 300)
		use_LU <- 1
	
	## Cholesky factorization
	
	L <- NULL
	resnrm <- NULL
	xx <- ones(m,1)
	if (!use_LU){
		nnzmat <- nnzero(coeff$mat11)
		nnzmatdiff <- (nnzmat != nnzmatold)   
		solve_ok <- 1
		solvesys <- 1 
		if((nnzmat > spdensity*m^2) | (m < 500)){  
			matfct_options <- "chol"
			if(is(schur,"sparseMatrix"))
				schur <- as.matrix(schur)
		}else{
			matfct_options <- "spchol"
			if(!is(schur,"sparseMatrix"))
				schur <- as(schur,"sparseMatrix")
		}
		if(printlevel)
			print(matfct_options)
		if(matfct_options=="chol"){
			L$matfct_options <- "chol"    
			L$perm <- c(1:m)
			L$L <- chol(schur) 
			
			if(class(L$L)=="try-error"){
				solve_ok <- -2
				solvesys <- 0
				print("  chol: Schur complement matrix not pos. def.")
			}
		}else if(matfct_options=="spchol"){
			L$matfct_options <- "spchol"    
			CX <- Cholesky(schur)
			if(class(CX)=="try-error"){
				solve_ok <- -2
				solvesys <- 0
				print("chol: Schur complement matrix not pos. def.")
			}
			r <- expand(CX)
			L$Lt <- r$L 
			L$L <- tt(r$L)
			L$perm <- r$P@perm
		}
		if(solvesys){
			if (ncolU){
				tmp <- tt(coeff$mat12)%*%linsysolvefun(L,coeff$mat12)-coeff$mat22 
				if(is(tmp,"sparseMatrix"))
					tmp <- as(tmp,"denseMatrix")
				elu <- expand(lu(tmp))
				L$Ml <- elu$L
				L$Mu <- elu$U
				L$Mp <- elu$P
				tol <- 1e-15 
				if(ncol(tmp)>1){
					max.diag <- max(abs(diag(as.matrix(L$Mu))))
					idx <- which(abs(diag(as.matrix(L$Mu)))/max.diag < tol)
					if(!is.null(idx)){
						pertdiag <- rep(0,ncolU)
						pertdiag[idx] <- tol 
						L$Mu <- L$Mu + max.diag*as(diag(pertdiag),"sparseMatrix")
					}
				}else{
					if(abs(L$Mu[1,1])<tol)
						L$Mu[1,1] = L$Mu[1,1] + tol
				}
			}
			qmr <- symqmr(coeff,rhs,L)
			xx <- qmr$xx
			resnrm <- qmr$resnrm
			solve_ok <- qmr$solve_ok
			if ((solve_ok<=0) & (printlevel))
				print(c("warning: symqmr fails:",solve_ok)) 
			if(printlevel>=3)
				print(length(resnrm)-1) 
		}
		if((solve_ok < 0)|(solvesys == 0)){
			if (m < 5000){ 
				use_LU <- 1
				if (printlevel)
					print("switch to LU factor") 
			}else if(solve_ok==-0.5){
				solve_ok <- 1 
			}
		}
	}
	## symmetric indefinite or LU factorization
	
	if(use_LU){
		nnzmat <- nnzero(coeff$mat11)+nnzero(coeff$mat12) 
		nnzmatdiff <- (nnzmat != nnzmatold)  
		solve_ok <- 1 
		if(!is.null(coeff$mat22))
			raugmat <- rBind(cBind(coeff$mat11,coeff$mat12),
					cBind(tt(coeff$mat12),coeff$mat22)) 
		else
			raugmat <- coeff$mat11 
		if((nnzmat > spdensity*m^2) | ((m+ncolU) < 500)) 
			matfct_options <- "lu"     
		else
			matfct_options <- "splu"
		if(printlevel)
			print(matfct_options)
		if(matfct_options=="lu"){
			if(is(raugmat,"sparseMatrix"))
				raugmat <- as(raugmat,"denseMatrix")
			lua = expand(lu(raugmat))
			L$l <- lua$L
			L$u <- lua$U
			L$p <- lua$P
			L$matfct_options <- "lu" 
			L$p = as(L$p,"sparseMatrix") 
			idx = which(abs(diag(L$u)) < 1e-15) 
			if(!is.null(idx)){
				if(printlevel)
					print("lu: modify diag(L$u)")
				n <- length(raugmat)
				lam <- zeros(n,1)
				lam[idx] <- 1e-15 
				L$u <- L$u + as(diag(lam),"sparseMatrix")
			}
		}
		if(matfct_options=="splu"){ 
			if(!is(raugmat,"sparseMatrix"))
				raugmat <- as(raugmat,"sparseMatrix")
			L$perm <- c(1:dim(raugmat)[1])
			L$matfct_options <- "splu"
			xA <- expand(lu(raugmat))
			L$l <- xA$L
			L$u <- xA$U
			L$p <- xA$P
			L$q <- xA$Q
		}
		xxA <- mybicgstab(coeff,rhs,L)
		xx <- xxA$xx
		resnrm <- xxA$resnrm
		solve_ok <- xxA$solve_ok
		if((solve_ok<=0) & (printlevel))
			print(c("warning: mybicgstab fails: ",solve_ok))
		if (printlevel>=3)
			print(length(resnrm)-1)
	}
	nnzmatold <- nnzmat 
	matfct_options_old <- matfct_options

        global_var$use_LU = use_LU
        global_var$solve_ok = solve_ok
        global_var$nnzmatold = nnzmatold
        global_var$matfct_options = matfct_options
        global_var$matfct_options_old = matfct_options_old

	return(list(xx=xx,coeff=coeff,L=L,global_var=global_var))
}
