sepelimdwd_weight = function(Xp,Xn,penalty,wgt){
	
	##Xp and Xn are matrices.  penalty is a scalar.  This is an
	##adaptation of the Matlab function of the same name, written
	##by J. S. Marron, available at
	##https://genome.unc.edu/pubsup/dwd/.
	
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
	XpnY11 <- XpnY[1,1]
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
	ym = y[2:n,1]
	
	## nv is the number of variables (eliminating beta)
	## nc is the number of constraints
	nv = 1 + dnew + 4*n
	nc = 2*n
	##Set up the block structure, constraint matrix, rhs, and cost vector
	
	blk = list()
	blk$type = character()
	blk$size = list()
	blk$type[1] = 'q'
	blk$size[[1]] = cbind(dnew+1,3*ones(1,n))
	blk$type[2] = 'l'
	blk$size[[2]] = n	
	
	Avec = list()
	A = zeros(nc,nv-n)
	col1 = RpnY[,1]
	A[1:(n-1),2:(dnew+1)] = t(RpnY[,2:n] - col1%*%t(ym))
	A[1:(n-1),seq(dnew+5,dnew+1+3*n,3)]	= -1*speye(n-1)
	A[1:(n-1),seq(dnew+6,dnew+2+3*n,3)]	= speye(n-1)
	A[1:(n-1),dnew+2] = ym
	A[1:(n-1),dnew+3] = -1*ym
	A[n,1] = 1
	A[(n+1):(n+n),seq(dnew+4,dnew+3+3*n,3)] = speye(n)
	
	Avec[[1]] = t(A)
	Avec[[2]] = (rbind(cbind(-1*ym,speye(n-1)),zeros(1+n,n)))
	b = rbind(zeros(n-1,1),ones(1+n,1))	
	if(is.null(wgt)){
		weight = rbind(2*np/n*ones(np,1),2*nn/n*ones(nn,1))
		weight = rbind(ones(np,1),ones(nn,1))
	}else{
		weight = rbind(wgt[1]*ones(np,1),wgt[2]*ones(nn,1))
        }
	C = list()
	c = zeros(nv-n,1)
	c[seq(dnew+2,dnew+1+3*n,3),1] = weight
	c[seq(dnew+3,dnew+2+3*n,3),1] = weight
		
	C[[1]] = c
	C[[2]] = penalty*weight
	
	##Solve the SOCP problem
	
	OPTIONS <- sqlparameters()
	spdensity <- NULL
	initial <- infeaspt(blk,Avec,C,b,spdensity=spdensity)
	X0 <- initial$X0
	lambda0 <- initial$y0
	Z0 <- initial$Z0
	##################################################################################
	## save data as example
	##################################################################################
#	nargin <<- 8
#	sqlpData = list(blk=blk,At=Avec,C=C,b=b,X0=X0,y0=lambda0,Z0=Z0)
#	save(sqlpData,file="DWDnew/data/sqlpData.rda")
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
	
	## Compute the normal vector w and constant term beta.
	
	barw = X1[2:(dnew+1)]
	if (d>n){
		w = Q %*% barw
	}else{
		w = barw	
	}
	beta = X1[dnew + 2] - X1[dnew + 3] - X2[1] - t(col1)%*%barw
	
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
	
	## Compute the minimum of the supposedly positive 
	## and the maximum of the supposedly negative residuals.
	## Refine the primal solution and print its objective value.
	
	residp = t(Xp) %*% w + beta[1] #optimization
	residn = t(Xn) %*% w + beta[1] #optimization
	minresidp = min(residp)
	maxresidn = max(residn)
	res = t(XpnY) %*% w + beta[1] * y
	rsc = 1/sqrt(penalty)
	xi = -1* res + rsc[1]
	xi[xi<0] <- 0
	
	totalviolation = sum(xi)
	minresidpmod = min(residp + xi[1:np])
	maxresidnmod = max(residn - xi[(np+1):n])
	minxi = min(xi)
	maxxi = max(xi)
	resn = res + xi
	rresn = 1 / resn
	primalobj = penalty * sum(weight*xi) + sum(weight*rresn)
	##print(primalobj)
	
	
	##Compute the dual solution alp and print its objective value.
	alp = zeros(n,1)
	lambda1 = lambda[1:(n-1)]
	alp[1] = -1*t(ym)%*%lambda1
	alp[2:n] = lambda1
	
	alp = alp * (as.numeric(alp>0)) + 1e-10
	
	sump = sum(alp[1:np])
	sumn = sum(alp[(np+1):n])
	sum2 = (sump + sumn)/2	
	alp[1:np] = (sum2/sump)*alp[1:np]
	alp[(np+1):n] = (sum2/sumn)*alp[(np+1):n]
	maxalp = max(alp/weight)
	if (maxalp > penalty | maxxi > 1e-3){
		alp = (penalty[1]/maxalp)*alp
	}
	minalp = min(alp)
	p = RpnY%*%alp
	
	eta = -1*normsvd(p)
	
	gamma = 2*sqrt(weight*alp)
	dualobj = eta + sum(gamma)
	
	##dualgap is a measure of the accuracy of the solution
	dualgap = abs(primalobj - dualobj)
	
	if (dualgap > 1e-4){
		flag = -1
		print(c(primalobj,dualobj))
	}
	
	return(list(w=w,beta=beta,obj=obj,residp=residp,residn=residn,
					alp=alp,totalviolation=totalviolation,
					dualgap=dualgap,flag=flag))
	
}

DWD1SM_weight = function(trainp,trainn,threshfact=100,wgt=NULL){
	
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
	sepelimout = sepelimdwd_weight(trainp,trainn,penalty,wgt)
	w = sepelimout$w
	beta = sepelimout$beta
	obj = sepelimout$obj
	flag = sepelimout$flag
	if (flag == -1){	
		cat("Inaccurate solution!\n")
	}
	if (flag == -2){
		cat("Infeasible or unbounded optimization problem!\n")
	}
	dirvec = w/normsvd(w)
	beta = beta/normsvd(w)
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

