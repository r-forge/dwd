###############################################################################
# TODO: DWD using unrestricted variables
# 
# Author: Hanwen Huang
###############################################################################

sepelimdwdHH = function(Xp,Xn,penalty){
	
	flag <- 0
	##Dimensions of the data
	dp <- dim(Xp)[1]
	np <- dim(Xp)[2]
	dn <- dim(Xn)[1]
	nn <- dim(Xn)[2]
	if (dn != dp) {stop('The dimensions are incomapatible.')}
	d <- dp
	
	##Dimension reduction in HDLSS setting.
	XpnY <- as.matrix(cbind(Xp,-1*Xn))
	n <- np + nn
	if(d>n){
		qrfact = qr(XpnY)
		Q = qr.Q(qrfact)
		R = qr.R(qrfact)
		RpnY = R
		dnew = n
	}else{
		RpnY = XpnY
		dnew = d
	}
	
	y = ones(np + nn,1)
	y[(np+1):(np+nn),1] = -1
	
	## nv is the number of variables (eliminating beta)
	## nc is the number of constraints
	nv = 2 + dnew + 4*n
	nc = 2*n + 1
	##Set up the block structure, constraint matrix, rhs, and cost vector
	
	blk = list()
	blk$type = character()
	blk$size = list()
	blk$type[1] = 'q'
	blk$size[[1]] = cbind(dnew+1,3*ones(1,n))
	blk$type[2] = 'l'
	blk$size[[2]] = n
	blk$type[[3]] = 'u'
	blk$size[[3]] = 1
	
	Avec = list()
	Aq = zeros(nc,1+dnew+3*n)
	Aq[1:n,2:(dnew+1)] = t(RpnY)
	
	Aq[1:n,seq(dnew+2,dnew+1+3*n,3)]	= -1*speye(n)
	Aq[1:n,seq(dnew+3,dnew+1+3*n,3)]	= speye(n)
	Aq[n+1,1] = 1
	Aq[(n+2):(n+n+1),seq(dnew+4,dnew+1+3*n,3)] = speye(n)
	Avec[[1]] = t(Aq)
	
	Al = rbind(speye(n),zeros(n+1,n))
	Avec[[2]] = Al
	
	Au = rbind(y,zeros(n+1,1))
	Avec[[3]] = Au
	
	b = rbind(zeros(n,1),ones(1+n,1))	
	
	C = list()
	c = zeros(1+dnew+3*n,1)
	c[seq(dnew+2,dnew+1+3*n,3),1] = ones(n,1)
	c[seq(dnew+3,dnew+1+3*n,3),1] = ones(n,1)
	
	
	C[[1]] = c
	C[[2]] = penalty*ones(n,1)
	C[[3]] = 0
	
	##Solve the SOCP problem
	
	library(Matrix)
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
	X3 <- X[[3]]    
	#####################################################
	
	## Compute the normal vector w and constant term beta.
	
	barw = X1[2:(dnew+1)]
	if (d>n){
		w = Q %*% barw
	}else{
		w = barw	
	}
	beta = X3
	
	normw = normsvd(w)
	if (normw < 1 - 1e-3){
		print(normw)
	}
	normwm1 = 0
	if (normw > 1 - 1e-3){
		w = w/normw
		normwm1 = normsvd(w) - 1
		beta = beta/normw
	}
	
	return(list(w=w,beta=beta,flag=flag,alp=lambda[1:n]))
	
}

dwdHH = function(trainp,trainn,threshfact=100){
	
	np = dim(trainp)[2]
	nn = dim(trainn)[2]
#	vpwdist2 = numeric(np*nn)
	vpwdist2x <- rdist(t(trainp),t(trainn))
#	for (ip in 1:np){
#		vpwdist2[((ip-1)*nn+1):(ip*nn)] <- colSums((trainp[,ip] - trainn)^2) #optimization
#	}
#	medianpwdist2 = median(vpwdist2)
	medianpwdist2 = median(vpwdist2x)^2
	
	penalty = threshfact / medianpwdist2
	sepelimout = sepelimdwdHH(trainp,trainn,penalty)
	w = sepelimout$w
	beta = sepelimout$beta
	flag = sepelimout$flag
	if (flag == -1){	
		cat("Inaccurate solution!\n")
	}
	if (flag == -2){
		cat("Infeasible or unbounded optimization problem!\n")
	}
	dirvec = w/normsvd(w)
	return(list(w=dirvec,beta=beta,alp=sepelimout$alp))
}

normsvd = function(aMatrix){
	##Reqturns the largest singular value of aMatrix
	o = svd(aMatrix,nu=0,nv=0)
	return(o$d[1])
}




