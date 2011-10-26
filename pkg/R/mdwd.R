###############################################################################
# TODO: Multiclass DWD
# 
# Author: Hanwen Huang
###############################################################################


sepelimmdwd = function(X,y,penalty,wgt){
	
	flag <- 0
	##Find the dimensions of the data
	d <- nrow(X)
	n <- ncol(X)
	nn = length(y)
	if (nn != n) {
		stop('The dimensions are incomapatible.')
	}
	K = length(unique(y))
	if(K<2)
		stop('Few than 2 classes.')			
	
	##Do the dimension reduction if in HDLSS setting.
	if(d>n){
		qrfact = qr(X)
		Q = qr.Q(qrfact)
		R = qr.R(qrfact)
		XX = R
		dnew = n
	}else{
		XX = X
		dnew = d
	}
	## nresid is the number of residuals,
	## nv is the number of variables and
	## nc is the number of constraints
	nresid = n*(K-1)
	nv = K + K*dnew + 1 + 3*nresid + nresid
	nc = 1 + 2*nresid + dnew + 1
	e = ones(nresid,1)
	speyen = speye(nresid)
	##Set up the block structure, constraint matrix, rhs, and cost vector.
	
	blk = list()
	blk$type = character()
	blk$size = list()
	blk$type[1] = 'q'
	blk$size[[1]] = cbind(K*dnew+1,3*ones(1,nresid))
	blk$type[2] = 'l'
	blk$size[[2]] = nresid
	blk$type[[3]] = 'u'
	blk$size[[3]] = K
	
	Avec = list()
	A1 = zeros(nc,1+K*dnew+3*nresid)
	A3 = zeros(nc,K)
	rind = 0
	rin2 = nresid+1
	cind = 1 + K*dnew
	for(k in 1:K){
		indk = which(y == k)
		if(is.null(indk)) 
			stop('Empty class.')
		Xk = XX[,indk]
		numk = length(indk)
		for(j in 1:K){
			if (j != k){
				A1[(rind+1):(rind+numk),((j-1)*dnew+2):(j*dnew+1)] = -t(Xk)
				A1[(rind+1):(rind+numk),((k-1)*dnew+2):(k*dnew+1)] =  t(Xk)
				A1[(rind+1):(rind+numk),seq(cind+1,cind+0+3*numk,by=3)] = -speye(numk)
				A1[(rind+1):(rind+numk),seq(cind+2,cind+1+3*numk,by=3)] =  speye(numk)
				A1[(rin2+1):(rin2+numk),seq(cind+3,cind+2+3*numk,by=3)] =  speye(numk)
				A3[(rind+1):(rind+numk),j] = -ones(numk,1)
				A3[(rind+1):(rind+numk),k] =  ones(numk,1)
				rind = rind + numk
				rin2 = rin2 + numk
				cind = cind + 3*numk
			}
		}
		A1[(1+2*nresid+2):nc,((k-1)*dnew+2):(k*dnew+1)] = speye(dnew)
	}
	A1[nresid+1,1] = 1
	A3[(1+2*nresid+1),] = ones(1,K)
	Avec[[1]] = t(A1)
	Avec[[2]] = t(rbind(speyen,zeros(1+nresid,nresid),zeros(1+dnew,nresid)))
	Avec[[3]] = t(A3)
	
	b = rbind(zeros(nresid,1),1,e,0,zeros(dnew,1))
	weighted_index = 1
	C = list()
	meann = n/K
	if(weighted_index==0){
		c = zeros(nv-nresid-K,1)
		c[seq(K*dnew+2,K*dnew+1+3*nresid,by=3)] = e
		c[seq(K*dnew+3,K*dnew+2+3*nresid,by=3)] = e
		C[[1]] = c
		C[[2]] = penalty*e
		C[[3]] = zeros(K,1)
	}else{
		weight = ones(K,1)
		c = ones(nv-nresid-K,1)
		c1 = ones(nresid,1)
		rind = K*dnew+1
		cind = 0
		for(k in 1:K){
			indk = which(y == k)
			if(is.null(indk)) stop('Empty class.')
			numk = length(indk)
			if(is.null(wgt))
				weight[k] = meann/numk
			else
				weight[k] = wgt[k]
			c[seq(rind+1,rind-2+3*(K-1)*numk,by=3)] = weight[k]*ones((K-1)*numk,1)
			c[seq(rind+2,rind-1+3*(K-1)*numk,by=3)] = weight[k]*ones((K-1)*numk,1)
			c1[(cind+1):(cind+(K-1)*numk)] = weight[k]*penalty*ones((K-1)*numk,1)
			rind = rind + 3*(K-1)*numk
			cind = cind + (K-1)*numk
		}
		C[[1]] = c
		C[[2]] = c1
		C[[3]] = zeros(K,1)
	}
	
	##Solve the SOCP problem
	
#	library(Matrix)
	spdensity <- NULL
	OPTIONS = sqlparameters()
	OPTIONS$maxit = 40
	OPTIONS$vers = 2
	OPTIONS$steptol = 1e-10
	OPTIONS$gaptol = 1e-12
	initial = infeaspt(blk,Avec,C,b,spdensity=spdensity)
	X0 = initial$X0
	lambda0 = initial$y0
	Z0 = initial$Z0
	soln = sqlp(blk,Avec,C,b,OPTIONS,X0,lambda0,Z0)
	obj = soln$obj
	X = soln$X
	lambda = soln$y
	Z = soln$Z
	info = soln$info
	##If infeasible or unbounded, break.
	
	if (info$termcode > 0) {
		flag = -2 
		return
	}
	
	##Compute the normal vectors W and constant terms beta.
	X1 = X[[1]]
	X2 = X[[2]] 
	beta = X[[3]]
	barW = matrix(X1[2:(K*dnew+1)],dnew,K)
	if (d>n)
		W = Q%*%barW
	else
		W = barW
	normw = normsvd(as.vector(W))
	if(normw < 1 - 1e-3) 
		print(normw) 
	normwm1 = 0
	if(normw > 1 - 1e-3){
		W = W / normw
		normwm1 = normsvd(as.vector(W))-1
		beta = beta / normw
	}
	
	## Compute the residuals.
	## Refine the primal solution and print its objective value.
	
	res = NULL
	for(k in 1:K){
		indk = which(y == k)
		if (length(indk) == 0) stop('Empty class.')
		Xk = XX[,indk]
		numk = length(indk)
		for(j in 1:K){
			if (j != k){
				res = rbind(res, t(Xk)%*%(barW[,k]-barW[,j]) + ones(numk,1)*(beta[k]-beta[j]))
			}
		}
	}
	rsc = 1/sqrt(penalty)
	xi = rsc - res
	xi[which(xi<0)] = 0
	totalviolation = sum(xi)
	minresidmod = min(res+xi)
	minxi = min(xi)
	maxxi = max(xi)
	resn = res + xi
	rresn = 1/resn
	primalobj = penalty*totalviolation + sum(rresn)
	if(weighted_index==1){
		ynew = NULL
		for(k in 1:K){
			indk = which(y==k)
			if (length(indk) == 0) stop('Empty class.')
			numk = length(indk)
			ynew = rbind(ynew,k*ones((K-1)*numk,1))
		}
		tempobj = penalty*xi + rresn
		primalobj = sum(tempobj*weight[ynew])
	}
	## Compute the dual solution alp and print its objective value.
	
	alp = lambda[1:nresid]
	alp[which(alp<0)] = 0
	
	maxalp = max(alp/weight[ynew])
	if (maxalp > penalty | maxxi > 1e-3)
		alp = (penalty/maxalp) * alp
	
	minalp = min(alp)
	p = t(A1[1:nresid,2:(K*dnew+1)])%*%alp
	eta = - normsvd(p)
	gamma = 2 * sqrt(alp*weight[ynew])
	dualobj = eta + sum(gamma)
	
	## dualgap is the duality gap, a measure of the accuracy of the solution.
	
	dualgap = primalobj - dualobj
	
	if (abs(dualgap) > 1e-3){
		flag = -1
		print(dualgap)
	}
	return(list(W=W,beta=beta,flag=flag))
	
}

mdwd = function(xMat,yLabel,threshfact=100,wgt=NULL){
	
	n = ncol(xMat)
	unilabel = unique(yLabel)
	K = length(unilabel)
	vpwdist = NULL
	for(i in 1:(K-1)){
		for(j in (i+1):K){
			dist <- rdist(t(xMat[,yLabel==i]),t(xMat[,yLabel==j]))
			vpwdist = c(vpwdist,as.vector(dist))
		}
	}
	medianpwdist2 = median(vpwdist^2)
	
	penalty = threshfact / medianpwdist2
	sepelimout = sepelimmdwd(xMat,yLabel,penalty,wgt)
	W = sepelimout$W
	beta = sepelimout$beta
	obj = sepelimout$obj
	flag = sepelimout$flag
	if (flag == -1){	
		cat("Inaccurate solution!\n")
	}
	if (flag == -2){
		cat("Infeasible or unbounded optimization problem!\n")
	}
	dirvec = W/normsvd(as.vector(W))
	beta = beta/normsvd(as.vector(W))
	return(list(w=dirvec,beta=beta,obj=obj))
}

normsvd = function(aMatrix){
	##Reqturns the largest singular value of aMatrix
	o = svd(aMatrix,nu=0,nv=0)
	return(o$d[1])
}

rdist = function(x1,x2){
	apply(x1,1,function(a){
				apply(x2,1,function(b){
							return(d=sqrt(sum((a-b)^2)))
						})
			})
}



