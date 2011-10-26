##*************************************************************************
## symqmr: symmetric QMR with left (symmetric) preconditioner. 
##         The preconditioner used is based on the analytical
##         expression of inv(A).  
##
## [x,resnrm,solve_ok] = symqmr(A,b,L,tol,maxit) 
##
## child function: linsysolvefun.m 
##
## A = [mat11 mat12; mat12' mat22].
## b = rhs vector.
## if matfct_options = 'chol' or 'spchol' 
##    L = Cholesky factorization of (1,1) block. 
##    M = Cholesky factorization of 
##        Schur complement of A ( = mat12'*inv(mat11)*mat12-mat22).
## else
##    L = triangular factors of A.
##    M = not relevant.
## end
## resnrm = norm of qmr-generated residual vector b-Ax. 
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##*************************************************************************

symqmr = function(A,b,L,tol,maxit) {
	
	resnrm <- NULL
	N <- length(b) 
	if (nargs() < 5)
		if(is.null(A$mat22)){
			maxit <- 5
		}else{
			maxit <- max(50,max(5,max(dim(A$mat22))))
		}
	if (nargs() < 4)
		tol <- 1e-10
	tolb <- min(1e-4,tol*normsvd(b))
	
	solve_ok <- 1 
	x <- zeros(N,1)
	Aq <- matvec(A,x)
	r <- b - Aq  
	err <- normsvd(r)
	resnrm[1] <- err
	minres <- err
	xx <- x 
	if (err < tolb)
		return
	q <- precond(A,L,r)
	tau_old <- normsvd(q)      
	rho_old <- t(r)%*%q 
	theta_old <- 0 
	d <- zeros(N,1) 
	res <- r
	Ad <- zeros(N,1)
	
	## main loop
	
	tiny <- 1e-30
	for(iter in 1:maxit){ 
		Aq <- matvec(A,q)
		sigma <- t(q)%*%Aq 
		if (abs(sigma) < tiny){
			solve_ok <- 2
			break
		}else{
			alpha <- rho_old/sigma 
			r <- r - alpha*Aq
		}
		u <- precond(A,L,r)
		theta <- normsvd(u)/tau_old
		c <- 1/sqrt(1+theta^2) 
		tau <- tau_old*theta*c
		gam <- (c^2*theta_old^2)
		eta <- (c^2*alpha) 
		d <- gam*d + eta[1,1]*q
		x <- x + d 
		
		Ad <- gam*Ad + eta[1,1]*Aq
		res <- res - Ad
		err <- normsvd(res)
		resnrm[iter+1] <- err 
		if (err < minres){
			xx <- x
			minres <- err
		}
		if (err < tolb) break
		if (iter > 5) {
			if (err > 0.98*mean(resnrm[(iter-5):iter])){
				solve_ok <- -0.5
				break 
			}
		}
		
		if (abs(rho_old) < tiny){
			solve_ok <- 2
			break
		}else{
			rho <- t(r)%*%u 
			beta <- rho/rho_old 
			q <- u + beta[1,1]*q 
		}
		rho_old <- rho 
		tau_old <- tau 
		theta_old <- theta
	}
	return(list(xx=xx,resnrm=resnrm,solve_ok=solve_ok))
}

## matvec: matrix-vector multiply.
## matrix = [A.mat11 A.mat12; t(A.mat12) A.mat22]

matvec = function(A,x){
	m <- max(dim(A$mat11))
	m2 <- length(x)-m 
	if (m2 > 0)
		x1 <- x[1:m] 
	else
		x1 <- x 
	Ax <- mexMatvec(A$mat11,x1)
	if (m2 > 0){
		x2 <- x[m+c(1:m2)]
		Ax <- Ax + mexMatvec(A$mat12,x2) 
		Ax2 <- mexMatvec(A$mat12,x1,1) + mexMatvec(A$mat22,x2)
		Ax <- c(Ax,Ax2)
	}
	return(Ax)
}

## precond: 

precond = function(A,L,x){
	m <- length(L$perm)
	m2 <- length(x)-m
	if (m2 > 0)
		x1 <- x[1:m] 
	else
		x1 <- x
	if (m2 > 0){
		x2 <- x[m+c(1:m2)]
		w <- linsysolvefun(L,x1) 
		z <- mexMatvec(A$mat12,w,1) - x2
		z <- solve(L$Mu)%*%solve(L$Ml)%*%z[L$Mp@perm]
		x1 <- x1 - mexMatvec(A$mat12,z) 
	}
	Mx <- linsysolvefun(L,x1)  
	if (m2 > 0)
		Mx = c(Mx,z[,1])
	return(Mx)
}

