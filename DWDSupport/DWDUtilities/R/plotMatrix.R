###############################################################################
## plotMatrix.R
##
## Function to plot projection matrix based on directions
## The number of directions can be either two or more
## The directions may not necessarily  orthogonal to each other
## 
## newest version for ploting graph matrix
##
## Author: Xiaosun Lu
################################################################################

## y value on one dimensional projection (like jittering)
height = function(x){
	num = length(x)
	d = max(density(x)$y) - min(density(x)$y)
	y1 = min(density(x)$y) + 0.3*d
	y2 = min(density(x)$y) + 0.7*d
	return( y1 + (1:num)/num * (y2-y1) )
}

vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)


plotMatrix2 = function(
		dirDf,
		featureDf=NULL,
		col ,## 
		shape = rep(1,nrow(featureDf)),
		alpha = 0.5,
		xlims = NULL){
	## xlims = a list of min, and max of x axis
	
	if( nrow(dirDf)!=ncol(featureDf) ){
		stop("Dimension does not match.")
	}
	if( length(col)!=nrow(featureDf) ){
		stop("Each sample should be assigned a color.")
	}
	if( length(shape)!=nrow(featureDf) ){
		stop("Each sample should be assigned a shape.")
	}
  #  browser()	
	
	labels = names(dirDf)
	#names(featureDf) 
	
	m = ncol(dirDf)
	#m
	# [1] 2
	
	if(length(alpha)==1){alpha=rep(alpha,nrow(featureDf))}
	########################################
	## project the raw data onto the directions
	projDf = data.frame(as.matrix(featureDf)%*%as.matrix(dirDf))
	head(projDf)
	#       OPC1    dwdDir
	# 1 6.171253 -1.278577
	# 2 5.678127 -1.386254
	# 3 5.670071 -1.199860
	# 4 5.543858 -1.381970
	# 5 6.170930 -1.228574
	# 6 6.687502 -1.550592
	
	
	projAllDf = projDf ## DF saved for density plot of mean values
	
	if(is.null(xlims)){
		dirMin = apply(projDf,2,min)
		dirMin
		#      OPC1    dwdDir 
		#  4.895814 -5.391162 
		
		dirMax = apply(projDf,2,max)
		temp = dirMax - dirMin
		dirMax = dirMax + temp/100
		dirMin = dirMin - temp/100
	}else{
		dirMin = xlims$min
		dirMax = xlims$max
	}
	
	########################################
	## plots
	grid.newpage()
	pushViewport(viewport(layout=grid.layout(m,m)))
	
	unortho = FALSE   ## ## there is no unorthogonal pair
	for(i in 1:m){   
		#################
		## diagonal plots
		## i = 1
		diri = labels[i]   ## direction i
		temp1Df = data.frame(x=projDf[,diri],y=height(projAllDf[,diri]))
		temp1Df = cbind(temp1Df,col,shape,alpha)
		dx = density(projAllDf[,diri])$x
		dy = density(projAllDf[,diri])$y
		temp2Df = data.frame(dx,dy)[dx>dirMin[i] & dx<dirMax[i],]
		g = ggplot(temp1Df,aes(x=x,y=y)) + geom_line(data=temp2Df,aes(dx,dy))+ 
				geom_point(aes(x=x,y=y,colour=I(col),shape=shape,alpha=alpha)) +
				scale_colour_identity()+
				opts(legend.position = "none")+ xlab(diri)+ ylab("density")+ xlim(dirMin[i],dirMax[i]) 
		print(g, vp=vplayout(i,i))
		
		##################
		## off-diagonal plots
		for(j in (1:m)[-i]){ 
			## j = 2
			diri = labels[i]
			diri
			# [1] "PC1"
			dirj = labels[j]
			dirj
			# [1] "PC2"
			
			############
			## the two directions are orthogonal 
			if(abs(t(dirDf[,i])%*%dirDf[,j])<1e-5 & j>=i+1){
				tempDf = data.frame(x=projDf[,dirj],y=projDf[,diri])
				tempDf = cbind(tempDf,col,shape,alpha)
				
				g = ggplot(tempDf,aes(x,y,colour=I(col),shape=shape)) + geom_point(aes(alpha=alpha)) + 
						scale_colour_identity()+
						opts(legend.position = "none")+ xlab(dirj)+ylab(diri)+
						xlim(dirMin[j],dirMax[j]) +ylim(dirMin[i],dirMax[i])
				print( g, vp=vplayout(i,j))
			}
			
			############
			## the two directions are not orthogonal
			if(abs(t(dirDf[,i])%*%dirDf[,j])>1e-5){
				unortho = TRUE   ## there is a unorthogonal pair
				d1 = dirDf[,i]  ## direction
				d2 = dirDf[,j]
				## ya is the vector orthogonal to xa
				## xa is the direction 1
				## b = -sum(d1^2)/(sum(d1*d2))
				b =(sum(d1*d2))
				xa = d1   ## x axis
				ya = d2-b*d1
				ya = ya/sqrt(sum(ya^2))
				
				## projection
				## d1 is one axis, but d2 is not orthogonal to d1, so x-value stays the same and y-value need to be adjusted
				x=c()
				y=c()
				xm = cbind(d1,d2)
				
				dmat = xm%*%solve(t(xm)%*%xm)%*%t(xm)
				vmat = dmat%*%t(featureDf)
				omat = cbind(xa,ya)
				xmat = t(t(omat) %*% vmat)
				
				## coordinates of projected points 
				tempDf = data.frame( x=xmat[,1],y=xmat[,2])
				xmin = min(tempDf[,"x"])
				xmax = max(tempDf[,"x"])
				ymin = min(tempDf[,"y"])
				ymax = max(tempDf[,"y"])
				
				x0 = sum(d2*xa)
				y0 = sum(d2*ya)
				
				tempDf = cbind(tempDf,col,shape,alpha)
				g = ggplot(tempDf,aes(x,y,colour=I(col),shape=shape)) + geom_point(aes(alpha=alpha)) + 
						scale_colour_identity()+
						geom_abline(intercept=0, slope=(y0/x0))+ geom_hline(yintercept=0)+
						opts(legend.position="none")+xlab(diri)+ylab(dirj)+
						xlim(dirMin[i],dirMax[i])
				print(g, vp=vplayout(j,i))		
			}
		}
	}
}
##=============================================================================
## color chart for plot matrix
##=============================================================================


colChart = function(color,class,savePlot=NULL,xlab){
	
	class = as.character(class)
	checkEquals(length(color),length(class))
	tab = table(class,color)
	temp = apply(tab,1,function(v){
				checkEquals(sum(v!=0),1)
				c(names(v)[which(v!=0)],v[which(v!=0)])
			})
	df = data.frame(n = 1:dim(tab)[1],class=colnames(temp),col=temp[1,],num = as.numeric(temp[2,]))
	rownames(df)=NULL
	df
	#   n   class       col num
	# 1 1  BBraun       red 375
	# 2 2  Baxter      blue 378
	# 3 3 Hospira darkgreen 405
	
	## plot
	if(!is.null(savePlot)){
		pdf(savePlot)
	}
	
	g = ggplot(df,aes(x=n, y=num)) + 
			geom_bar(stat="identity",aes(colour=I(col),fill=I(col)))+
			scale_colour_identity()+
			scale_fill_identity()+
			opts(legend.position="none")+
			xlab(xlab)+ylab("Number of Samples")+
			geom_text(aes(x=n, y=rep(max(num)/2,length(num)), 
							label=class),angle=0)+
			coord_flip()
	print(g)
	
	if(!is.null(savePlot)){
		dev.off()
	}
	
	
	return(df)
}









