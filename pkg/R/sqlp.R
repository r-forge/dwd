
## sqlp: solve an SOCP program by infeasible path-following method. 
##
##  [obj,X,y,Z,info,runhist] = ...
##       sqlp(blk,At,C,b,OPTIONS,X0,y0,Z0);
##
##  Input: blk: a cell array describing the structure of the SOCP data.
##          At: a cell array with At{p} 
##         b,C: data for the SOCP instance.
##  (X0,y0,Z0): an initial iterate (if it is not given, the default is used).
##     OPTIONS: a structure that specifies parameters required in sqlp.m,
##              (if it is not given, the default in sqlparameters.m is used). 
##
##  Output: obj  = [<C,X> <b,y>].
##          (X,y,Z): an approximately optimal solution or a primal or dual
##                   infeasibility certificate. 
##          info.termcode = termination-code  
##          info.iter     = number of iterations
##          info.cputime  = total-time
##          info.gap      = gap
##          info.pinfeas  = primal_infeas
##          info.dinfeas  = dual_infeas  
##          runhist.pobj    = history of primal objective value. 
##          runhist.dobj    = history of dual   objective value.
##          runhist.gap     = history of <X,Z>. 
##          runhist.pinfeas = history of primal infeasibility. 
##          runhist.dinfeas = history of dual   infeasibility. 
##          runhist.cputime = history of cputime spent.
##          (Xiter,yiter,Ziter): last iterates.
##
##*************************************************************************
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
#%*****************************************************************************

sqlp = function(blk,At,C,b,OPTIONS,X0,y0,Z0){
	
        global_var = NULL
	use_LU <- 0
	nnzmatold <- NULL
	
	vers        <- 1
	predcorr    <- 1
	gam         <- 0
	expon       <- 1
	gaptol      <- 1e-8
	inftol      <- 1e-8
	steptol     <- 1e-6
	maxit       <- 100
	printlevel  <- 0
	stoplevel   <- 1
	plotyes     <- 1
	spdensity   <- 0.5
	rmdepconstr <- 0
	cachesize   <- 256
	
	if(!missing(OPTIONS)){
		vers <- OPTIONS$vers
		predcorr <- OPTIONS$predcorr
		gam <- OPTIONS$gam
		expon <- OPTIONS$expon
		gaptol <- OPTIONS$gaptol
		inftol <- OPTIONS$inftol
		steptol <- OPTIONS$steptol
		maxit <- OPTIONS$maxit
		printlevel <- OPTIONS$printlevel
		plotyes <- OPTIONS$plotyes
		spdensity <- OPTIONS$spdensity
		rmdepconstr <- OPTIONS$rmdepconstr
		cachesize <- OPTIONS$cachesize
	}
	X <- X0
	y <- y0
	Z <- Z0  
	
        global_var$use_LU = use_LU
        global_var$printlevel = printlevel
        global_var$spdensity = spdensity
        global_var$nnzmatold = nnzmatold
        global_var$solve_ok = NULL
        global_var$matfct_options = NULL
        global_var$matfct_options_old = NULL

	# validate SOCP data. 
	
	tstart <- proc.time()[3] 
	validated <- validate(blk,At,C,b,X,y,Z,spdensity=spdensity)
	At <- validated$At
	C <- validated$C
	dim <- validated$dim
	numblk <- validated$numblk
	X <- validated$X
	Z <- validated$Z
        global_var$spdensity = validated$spdensity

	# convert unrestricted blk to
	# linear blk.
	
	ublkidx <- zeros(length(blk$type),1)
	for(p in 1:length(blk$type)){ 
		if(blk$type[p]=="u" & sum(blk$size[[p]]) > 20){
			ublkidx[p] <- 1 
			n <- 2*blk$size[[p]] 
			blk$type[p] <- "l" 
			blk$size[[p]] <- n
			At[[p]] <- rbind(At[[p]],-At[[p]])       
			C[[p]] <- rbind(C[[p]],-C[[p]])
			b2 <- 1 + abs(t(b))  
			normC <- 1+normsvd(C[[p]])
			normA <- 1+sqrt(apply(At[[p]]*At[[p]],2,sum))
			X[[p]] <- max(1,max(b2/normA))*ones(n,1)
			Z[[p]] <- max(1,max(c(normA,normC))/sqrt(n))*ones(n,1)
		}
	}
	
	# check if the matrices Ak are
	# linearly independent.
	
	m0 <- length(b)
	
	checklinear <- checkdepconstr(blk,At,b,y,rmdepconstr)
	At <- checklinear$At
	b <- checklinear$b
	y <- checklinear$y
	indeprows <- checklinear$idxB
	depconstr <- checklinear$ndepconstr
        global_var$depconstr = depconstr
	feasible <- checklinear$feasible
	if (!feasible){
		print("sqlp: SOCP is not feasible")
		return
	}
	
	# initialization
	A <- list() 
	for(p in 1:length(blk$type)){
		A[[p]] <- tt(At[[p]])
	}
	
	normb <- normsvd(b)
	
	normC <- ops(C,'norm')
	normA <- ops(A,'norm') 
	normX0 <- ops(X0,'norm')
	normZ0 <- ops(Z0,'norm') 
	m <- length(b)
	n <- ops(C,'getM') 
	AX <- AXfun(blk,A,X) 
	
	rp <- b-AX
	
	Aty <- Atyfun(blk,At,y)
	ZpATy <- ops(Z,'+',Aty)
#	nargin <<- 2
	ZpATynorm <- ops(ZpATy,'norm')
#	nargin <<- 3
	Rd <- ops(C,'-',ZpATy)
	obj <- c((blktrace(blk,C,X))[1,1],tt(b)%*%y)
	gap <- blktrace(blk,X,Z)[1,1]  
	mu <- gap/n  
	rel_gap <- gap/(1+mean(abs(obj)))
	prim_infeas <- normsvd(rp)/(1+normb)
#	nargin <<- 2
	dual_infeas <- ops(Rd,'norm')/(1+normC)
	infeas_meas <- max(prim_infeas,dual_infeas)
	pstep <- 0
	dstep <- 0
	pred_convg_rate <- 1
	corr_convg_rate <- 1
	prim_infeas_bad <- 0
	termcode <- -6
	runhist <- list()
	runhist$pobj <- obj[1]
	runhist$dobj <- obj[2] 
	runhist$gap <- gap
	runhist$pinfeas <- prim_infeas
	runhist$dinfeas <- dual_infeas
	runhist$infeas <- infeas_meas  
	runhist$step <- 0 
	cputime <- proc.time()[3] 
	runhist$cputime <- cputime-tstart 
	ttime <- list()
	ttime$preproc <- runhist$cputime 
	ttime$pred <- 0
	ttime$predstep <- 0 
	ttime$corr <- 0
	ttime$corrstep <- 0
	ttime$misc <- 0
	
	# display parameters, and
	# initial info
	
	if(printlevel>=2){
		print(c("num. of constraints = ",m))
		if(dim[1]){ 
			print(c(" dim. of socp   var  = ",dim[1])) 
			print(c(" num. of socp blk  = ",numblk[1])) 
		}
		if(dim[2])
			print(c(" dim. of linear var  = ",dim[2]))
		if(dim[3])
			print(c(" dim. of free   var  = ",dim[3]))
		print(" SDPT3: Infeasible path-following algorithms")
		if (printlevel>=3){       
			timed = mytimed(ttime$preproc)
			hh <- timed$hh
			mm <- timed$mm
			ss <- timed$ss
			print(" version  predcorr  gam  expon")
			print(c(" NT ",predcorr,gam,expon))
			print("  pstep dstep p_infeas d_infeas  gap mean(obj)  cputime")
			print(c(0,0,0,prim_infeas,dual_infeas,gap,mean(obj),hh,mm,ss))
		}
	}
	
	# start main loop
	
	indef <- rep(0,2)
	Xcholf <- blkcholfun(blk,X) 
	Xchol <- Xcholf$Xchol
	indef[1] <- Xcholf$indef
	Zcholf <- blkcholfun(blk,Z)
	Zchol <- Zcholf$Xchol
	indef[2] <- Zcholf$indef
	
	if(any(as.integer(indef))){
		if (printlevel)
			print("Stop: X or Z not positive definite")
		termcode <- -3
		return
	}
	mupredhist <- NULL
	for(iter in 1:maxit){  
                global_var$iter = iter
		
		update_iter <- 0
		breakyes <- 0
		pred_slow <- 0
		corr_slow <- 0
		step_short <- 0 
		tstart <- proc.time()[3]  
		time <- zeros(1,11) 
		time[1] <- proc.time()[3]
		
		# predictor step.
		
		if (predcorr)
			sigma <- 0 
		else{ 
			sigma <- 1-0.9*min(pstep,dstep) 
			if (iter == 1)
				sigma <- 0.5
		}
		sigmu <- sigma*mu
		ntpredict <- NTpred(blk,A,rp,Rd,sigmu,X,Z,global_var)
		par <- ntpredict$par
		dX <- ntpredict$dX
		dy <- ntpredict$dy
		dZ <- ntpredict$dZ
		coeff <- ntpredict$coeff
		L <- ntpredict$L
		hRd <- ntpredict$hRd
                global_var = ntpredict$global_var
                solve_ok = global_var$solve_ok
		
		if (solve_ok <= 0){
			runhist$cputime[iter+1] <- proc.time()[3]-tstart
			termcode <- -4
		}
		time[2] <- proc.time()[3]
		ttime$pred <- ttime$pred + time[2]-time[1]
		
		# step-lengths for predictor
		# step
		
		if (gam == 0) 
			gamused <- 0.9 + 0.09*min(pstep,dstep) 
		else
			gamused <- gam
		Xstep <- steplength(blk,X,dX) 
		if((Xstep>.99e12)&(blktrace(blk,C,dX)[1,1] < -1e-3)&(prim_infeas<1e-3)){
			if(printlevel)
				print(" Predictor: dual seems infeasible.")
		}
		pstep <- min(1,abs(gamused*Xstep))
		Zstep <- steplength(blk,Z,dZ) 
		time[3] <- proc.time()[3]        
		if((Zstep>.99e12)&(t(b)%*%dy>1e-3)&(dual_infeas<1e-3)){
			if(printlevel)
				print("Predictor: primal seems infeasible.")
		}
		dstep <- min(1,abs(gamused*Zstep))
		gappred <- blktrace(blk,ops(X,"+",dX,pstep),ops(Z,"+",dZ,dstep))
		mupred <- gappred/n
		mupredhist[iter] <- mupred[1,1] 
		ttime$predstep <- ttime$predstep + time[3]-time[2]
		
		#  stopping criteria for
		#  predictor step.
		
		if ((min(abs(pstep),abs(dstep))<steptol) & (stoplevel>0)){
			if (printlevel) {
				print("  Stop: steps in predictor too short:")
				print(c(" pstep = ",pstep,  "dstep = ",dstep))
			}
			runhist$cputime[iter+1] <- proc.time()[3]-tstart 
			termcode <- -2 
			breakyes <- 1 
		}
		
		if (iter >= 2){ 
			idx <- c(max(2,iter-2):iter)
			pred_slow <- all(mupredhist[idx]/mupredhist[idx-1]>0.4)
			idx <- c(max(2,iter-5):iter)
			pred_convg_rate <- mean(mupredhist[idx]/mupredhist[idx-1])
			pred_slow <- pred_slow + (mupred/mu > 5*pred_convg_rate)
		}
		
		if(!predcorr){
			if((max(mu,infeas_meas)<1e-6) & (pred_slow) & (stoplevel)){
				if (printlevel){ 
					print("lack of progress in predictor:")
					print(c("mupred/mu =",mupred/mu,
									"pred_convg_rate =",pred_convg_rate))
				}
				runhist$cputime[iter+1] <- proc.time()[3]-tstart 
				termcode <- -1 
				breakyes <- 1
			}else{ 
				update_iter <- 1 
			}
		}
		
		# corrector step.
		mupred <- mupred[1,1]
		if ((predcorr) & (!breakyes)){
			step_pred <- min(pstep,dstep)
			if (mu > 1e-6){
				if (step_pred < 1/sqrt(3)){ 
					expon_used <- 1
				}else{
					expon_used <- max(expon,3*step_pred^2) 
				}
			}else{ 
				expon_used <- max(1,min(expon,3*step_pred^2)) 
			}
			sigma <- min( 1, (mupred/mu)^expon_used ) 
			sigmu <- sigma*mu 
			
			ntcorrt <- NTcorr(blk,A,par,rp,Rd,sigmu,hRd,dX,dZ,coeff,L,X,Z,global_var)
			dX <- ntcorrt$dX
			dy <- ntcorrt$dy
			dZ <- ntcorrt$dZ
                        global_var = ntcorrt$global_var
                        solve_ok = global_var$solve_ok
			if (solve_ok <= 0){
				runhist$cputime[iter+1] <- proc.time()[3]-tstart 
				termcode <- -4
			}
			time[4] <- proc.time()[3]
			ttime$corr <- ttime$corr + time[4]-time[3]
			
			# step-lengths for corrector
			# step
			
			if (gam == 0) 
				gamused <- 0.9 + 0.09*min(pstep,dstep) 
			else
				gamused <- gam
			Xstep <- steplength(blk,X,dX)
			if((Xstep > .99e12) & (blktrace(blk,C,dX)[1,1] < -1e-3)
					& (prim_infeas < 1e-3)){
				pstep <- abs(Xstep)
				if (printlevel)
					print("Corrector: dual seems infeasible.")
			}else
				pstep <- min(1,abs(gamused*Xstep))
			Zstep <- steplength(blk,Z,dZ)
			time[5] <- proc.time()[3]
			if((Zstep>.99e12)&(t(b)%*%dy>1e-3)&(dual_infeas<1e-3)){
				dstep <- abs(Zstep)
				if (printlevel)
					print(" Corrector: primal seems infeasible.")
			}else
				dstep <- min(1,abs(gamused*Zstep))
#			nargin <<- 4
			gapcorr <- blktrace(blk,ops(X,'+',dX,pstep),ops(Z,'+',dZ,dstep)) 
			mucorr <- gapcorr/n
			ttime$corrstep <- ttime$corrstep + time[5]-time[4]
			
			
			#  stopping criteria for
			#  corrector step
			
			if (iter >= 2){ 
				idx <- c(max(2,iter-2):iter)
				corr_slow <- all(runhist$gap[idx]/runhist$gap[idx-1]>0.8) 
				idx <-c(max(2,iter-5):iter)
				corr_convg_rate <- mean(runhist$gap[idx]/runhist$gap[idx-1])
				corr_slow <- corr_slow+(mucorr/mu>max(min(1,5*corr_convg_rate),0.8)) 
			}
			if((max(mu,infeas_meas)<1e-6)&(iter>10)&(sum(corr_slow))&(stoplevel)){
				if(printlevel){ 
					print("lack of progress in corrector:")
					print(c("mucorr/mu =",mucorr/mu, "corr_convg_rate =",
									corr_convg_rate))
				}
				runhist$cputime[iter+1] <- proc.time()[3]-tstart 
				termcode <- -1 
				breakyes <- 1
			}else{
				update_iter <- 1
			}
		}  
		# udpate iterate
		
		indef <- c(1,1)
		if (update_iter){
			for(t in 1:10){
				Xcholf <- blkcholfun(blk,ops(X,'+',dX,pstep))
				Xchol <- Xcholf$chol
				indef[1] <- Xcholf$indef
				if (indef[1])
					pstep <- 0.8*pstep
				else
					break                        
			}
			for(t in 1:10){
				Zcholf <- blkcholfun(blk,ops(Z,'+',dZ,dstep))
				Zchol <- Zcholf$chol
				indef[2] <- Zcholf$indef
				if(indef[2])
					dstep <- 0.8*dstep
				else
					break                        
			}
			AXtmp <- AX + pstep*AXfun(blk,A,dX)
			prim_infeasnew <- normsvd(b-AXtmp)/(1+normb)
			if(any(as.integer(indef))){
				if (printlevel){
					print(" Stop: X, Z not both positive definite")
				}
				termcode <- -3
				breakyes <- 1         
			}else if((prim_infeasnew> max(c(rel_gap,50*prim_infeas,1e-8)))
					& (max(pstep,dstep)<=1) & (stoplevel)){
				if (printlevel)
					print(c("Stop: primal infeas deteriorated too much",
									prim_infeasnew))
				termcode <- -7 
				breakyes <- 1
			}else{
				X <- ops(X,'+',dX,pstep)  
				y <- y + dstep*dy           
				Z <- ops(Z,'+',dZ,dstep)  
			}
		}
		
		# adjust linear blk arising
		# from unrestricted blk
		
		for(p in 1:length(blk$type)){
			if (ublkidx[p]==1){
				len <- blk$size[[p]]/2
				alpha <- 0.8
				xtmp <- apply(cbind(X[[p]][1:len],X[[p]][len+c(1:len)]),2,min) 
				X[[p]][1:len] <- X[[p]][1:len] - alpha*xtmp
				X[[p]][len+c(1:len)] <- X[[p]][len+c(1:len)] - alpha*xtmp
				if (mu < 1e-4)
					Z[[p]] <- 0.5*mu/apply(cbind(X[[p]],rep(1,length(X[[p]]))),2)
				else{
					ztmp0 <- apply(cbind(Z[[p]][1:len],Z[[p]][len+c(1:len)]),2,max)
					ztmp <- apply(cbind(ztmp0,rep(1,length(ztmp0))),2,min) 
					beta1 <- t(xtmp)%*%(Z[[p]][1:len]+Z[[p]][len+c(1:len)])
					beta2 <- t(X[[p]][1:len]+X[[p]][len+c(1:len)]-2*xtmp)%*%ztmp
					beta <- max(0.1,min(beta1/beta2,0.5))
					Z[[p]][1:len] <- Z[[p]][1:len] + beta*ztmp
					Z[[p]][len+c(1:len)] <- Z[[p]][len+c(1:len)] + beta*ztmp
				}
			}
		}
		
		# compute rp, Rd,
		# infeasibities, etc.
		
		gap <- blktrace(blk,X,Z)[1,1]
		mu <- gap/n
		AX <- AXfun(blk,A,X) 
		rp <- b-AX
#		nargin <<- 3
		ZpATy <- ops(Z,'+',Atyfun(blk,At,y))
#		nargin <<- 2
		ZpATynorm <- ops(ZpATy,'norm')
#		nargin <<- 3
		Rd <- ops(C,'-',ZpATy)
		obj <- c(blktrace(blk,C,X)[1,1],  t(b)%*%y)
		rel_gap <- gap/(1+mean(abs(obj)))
		prim_infeas <- normsvd(rp)/(1+normb)
#		nargin <<- 2
		dual_infeas <- ops(Rd,'norm')/(1+normC)
		infeas_meas <- max(prim_infeas,dual_infeas)
		runhist$pobj[iter+1] <- obj[1]
		runhist$dobj[iter+1] <- obj[2]
		runhist$gap[iter+1] <- gap
		runhist$pinfeas[iter+1] <- prim_infeas
		runhist$dinfeas[iter+1] <- dual_infeas
		runhist$infeas[iter+1] <- infeas_meas
		runhist$step[iter+1] <- min(pstep,dstep) 
		runhist$cputime[iter+1] <- proc.time()[3]-tstart 
		time[6] <- proc.time()[3]
		ttime$misc <- ttime$misc + time[6]-time[5] 
		if(printlevel>=3){
			timed <- mytimed(sum(runhist$cputime))
			hh <- timed$hh
			mm <- timed$mm
			ss <- timed$ss
			print(c(iter,pstep,dstep))
			print(c(prim_infeas,dual_infeas,gap))
			print(c(mean(obj),hh,mm,ss))
		}
		
		# check convergence.
		
		if((ops(X,'norm')>1e15*normX0)|(ops(Z,'norm')>1e15*normZ0)){
			termcode <- 3
			breakyes <- 1 
		}
		if(obj[2] > ZpATynorm / max(inftol,1e-13)){
			termcode <- 1
			breakyes <- 1
		}
		if (-obj[1] > normsvd(AX) / max(inftol,1e-13)){
			termcode <- 2
			breakyes <- 1
		}
		if (max(rel_gap,infeas_meas) < gaptol){
			if (printlevel)
				print(c("Stop: max(relative gap, infeasibilities) <",gaptol))
			termcode <- 0
			breakyes <- 1
		}
		if(stoplevel){
			min_prim_infeas <- min(runhist$pinfeas[1:iter]) 
			prim_infeas_bad <- prim_infeas_bad + 
					(prim_infeas>max(1e-10,min_prim_infeas) & (min_prim_infeas<1e-2))
			if(mu < 1e-8)
				idx = c(max(1,iter-1):iter)
			else if(mu < 1e-4)
				idx = c(max(1,iter-2):iter) 
			else
				idx = c(max(1,iter-3):iter)
			idx2 <- c(max(1,iter-4):iter) 
			gap_ratio2 <- runhist$gap[idx2+1]/runhist$gap[idx2]
			gap_slowrate <- min(0.8,max(0.6,2*mean(gap_ratio2)))
			gap_ratio <- runhist$gap[idx+1]/runhist$gap[idx] 
			if((infeas_meas<1e-4|prim_infeas_bad)&(rel_gap<5e-3)){ 
				gap_slow <- all(gap_ratio>gap_slowrate) & (rel_gap<5e-3)
				if(vers==1) 
					tmptol <- max(prim_infeas,1e-2*dual_infeas) 
				else
					tmptol <- max(0.5*prim_infeas,1e-2*dual_infeas) 
				if (rel_gap < tmptol){ 
					if (printlevel)
						print("Stop: relative gap < infeasibility.")
					termcode <- 0
					breakyes <- 1           
				}else if (gap_slow){ 
					if (printlevel)
						print(" Stop: progress is too slow.")
					termcode <- -5 
					breakyes <- 1
				}
			}else if((prim_infeas_bad)&(iter>50)&all(gap_ratio>gap_slowrate)){
				if (printlevel)
					print(" Stop: progress is bad.")
				termcode <- -5
				breakyes <- 1 
			}else if((infeas_meas<1e-8)&(gap>1.2*mean(runhist$gap[idx]))){
				if (printlevel)
					print(" Stop: progress is bad.")
				termcode <- -5
				breakyes <- 1  
			}
			if((max(runhist$infeas)>1e-4)&
					(min(runhist$infeas)<1e-4|prim_infeas_bad)){ 
				rel_gap2 <- abs(diff(obj))/(1+mean(abs(obj))) 
				if(rel_gap2 < 1e-3) 
					step_short <- all(runhist$step[iter:(iter+1)]<0.1)
				else if(rel_gap2 < 1) {
					idx <- c(max(1,iter-3):(iter+1))
					step_short <- all(runhist$step[idx] < 0.05) 
				}
				if(step_short){ 
					if(printlevel)
						print(" Stop: steps too short consecutively")
					termcode <- -5
					breakyes <- 1      
				}
			}
		}
		if(breakyes)
			break
	}
	
	# end of main loop
	
	if ((termcode == -6) & (printlevel))
		print(" Stop: maximum number of iterations reached.")
	
	# produce infeasibility
	# certificates if appropriate
	
	if (iter >= 1){
		param <- list()
		param$obj         <- obj
		param$rel_gap     <- rel_gap 
		param$prim_infeas <- prim_infeas
		param$dual_infeas <- dual_infeas
		param$inftol      <- inftol
		param$m0          <- m0
		param$indeprows   <- indeprows
		param$termcode    <- termcode
		param$AX          <- AX 
		param$normX0      <- normX0 
		param$normZ0      <- normZ0 
		misc <- sqlpmisc(blk,A,At,C,b,X,y,Z,param)
		X <- misc$X
		y <- misc$y
		Z <- misc$Z
		resid <- misc$resid
		reldist <- misc$reldist
	}
	
	# recover unrestricted blk
	# from linear blk
	
	for(p in 1:length(blk$type)){
		if(ublkidx[p] == 1){
			n <- blk$size[[p]]/2 
			X[[p]] <- X[[p]][1:n]-X[[p]][n+c(1:n)] 
			Z[[p]] <- Z[[p]][1:n] 
		}
	}
	
	#  print summary
	info <- list()
	info$termcode <- termcode
	info$iter <- iter 
	info$gap <- gap 
	info$pinfeas <- prim_infeas
	info$dinfeas <- dual_infeas
	info$cputime <- sum(runhist$cputime) 
	if ((termcode == 1) | (termcode == 2))
		info$resid <- resid
	
	nnorm <- list()
	nnorm$b <- normb
	nnorm$C <- normC
	nnorm$A <- normA 
	nnorm$X <- ops(X,'norm')
	nnorm$y <- normsvd(y)
	nnorm$Z <- ops(Z,'norm') 
	sqlpsummary(runhist,ttime,termcode,resid,reldist,nnorm,global_var)
	return(list(obj=obj,X=X,y=y,Z=Z,info=info,runhist=runhist))  
}

