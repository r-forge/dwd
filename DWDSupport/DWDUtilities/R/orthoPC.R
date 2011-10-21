###############################################################################
## orthoPC.R
##
## Function to compute OPC's orthogonal to one known direction
## 
## Author: Xiaosun Lu
###############################################################################


orthoPC = function(dataMat, dir){
	
	## check the inputs
	if( !(length(dir)%in%dim(dataMat))){
		stop("Dimension does not match.")
	}
	if( !is.vector(dir) ){
		warning("dir should be a vector")
	}
	
	dir = dir/sqrt(sum(dir^2))
	## space orthogonal to DWD dir
	xm = dir
	ym = t(dataMat)
	## project ym onto xm
	projMat = xm%*%solve(t(xm)%*%xm)%*%t(xm)%*%ym
	projMat = t(projMat)
	dim(projMat)
	
	## projection on orthogonal space
	oprojMat = dataMat-projMat
	ofit = prcomp(oprojMat)
	dirDf = data.frame(dir,ofit$rotation)
	npc = min(ncol(dataMat),nrow(dataMat))
	names(dirDf)[-1] = paste("OPC",1:npc,sep="")
	head(dirDf)
	
	if(nrow(dataMat)>ncol(dataMat)){
		checkEquals(abs(sum(dirDf$dir%*%dirDf[,npc+1]))>0.9,TRUE)
		dirDf = dirDf[,-(npc+1)]
		head(dirDf)
		#                     dir        OPC1       OPC2       OPC3
		# Sepal.Length -0.2033417 0.706774863  0.6672503  0.1178917
		# Sepal.Width   0.2966902 0.703541803 -0.6410536 -0.0778081
		# Petal.Length -0.8638745 0.074094837 -0.2674790 -0.4203401
		# Petal.Width  -0.3526300 0.002859955 -0.2688535  0.8963044
	}
	
	return(dirDf)
}


