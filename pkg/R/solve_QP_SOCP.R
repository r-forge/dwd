solve_QP_SOCP <- function(Dmat,dvec,Amat,bvec){
	d <- ncol(Dmat)
	n <- ncol(Amat)
	blk <- list()
	blk$type <- character()
	blk$size <- list()
	blk$type[1] <- 'q'
	blk$size[[1]] <- d + 1
	blk$type[2] <- 'l'
	blk$size[[2]] <- n
	R <- chol(Dmat)
	Avec = list()
	A = matrix(0,n,d+1)
	A[,2:(d+1)] = t(Amat)%*%solve(R)
	Avec[[1]] = A
	Avec[[2]] <- -diag(n)
	b = bvec - t(Amat)%*%solve(Dmat)%*%dvec
	C = list()
	c = matrix(0,d+1,1)
	c[1] = 1
	C[[1]] = c
	C[[2]] = matrix(0,n,1)
	OPTIONS <- sqlparameters()
	spdensity <- NULL
	initial <- infeaspt(blk,Avec,C,b,spdensity=spdensity)
	X0 <- initial$X0
	lambda0 <- initial$y0
	Z0 <- initial$Z0
#	nargin <<- 8
	soln <- sqlp(blk,Avec,C,b,OPTIONS,X0,lambda0,Z0)
	obj <- soln$obj
	X <- soln$X
	lambda <- soln$y
	Z <- soln$Z
	info <- soln$info
	if(info$termcode>0){
		flag <- -2
		return
	}
	X1 <- X[[1]]
	X2 <- X[[2]]    
	w = X1[2:(d+1)]
	solution <- solve(R)%*%w + solve(Dmat)%*%dvec
	return(list(solution=solution))
}

