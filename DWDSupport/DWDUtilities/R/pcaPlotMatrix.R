###############################################################################
## pcaPlotMatrix.R
##
## pca plot matrix, which is good for showing the original data with dim > 2
## 
## Author: Xiaosun Lu
###############################################################################


pcaPlotMatrix = function(
		dataDf,
		classColors=NULL,
		npc=NULL,
		alpha=1){
	
	# head(dataDf)
	#   Sepal.Length Sepal.Width Petal.Length Petal.Width  class
	# 1          5.1         3.5          1.4         0.2 setosa
	# 2          4.9         3.0          1.4         0.2 setosa
	# 3          4.7         3.2          1.3         0.2 setosa
	# 4          4.6         3.1          1.5         0.2 setosa
	# 5          5.0         3.6          1.4         0.2 setosa
	# 6          5.4         3.9          1.7         0.4 setosa
	
	
	if( !("class"%in%names(dataDf)) ){
		stop("dataDf should have a class column.")
	}
	
	## PCA
	dataMat = as.matrix(dataDf[,names(dataDf)!="class"])
	fit = prcomp(dataMat)
	dirDf = data.frame(fit$rotation)
	if(is.null(npc)){ npc = min(ncol(dataMat),nrow(dataMat))}
	dirDf = dirDf[,1:npc]
	plotMatrix2(
			dirDf = dirDf,
			featureDf = scale(as.data.frame(dataMat),center=TRUE,scale=FALSE),
			col = classColors[as.numeric(dataDf$class)],
			alpha=alpha
	)	
	
}




##=============================================================================
## original data plots
##=============================================================================

originPlotMatrix = function(
		dataDf,
		#classes=NULL,
		classColors=NULL,
		nvars=NULL,
		alpha=1){
	
#	browser()

	# head(dataDf)
	#   Sepal.Length Sepal.Width Petal.Length Petal.Width  class
	# 1          5.1         3.5          1.4         0.2 setosa
	# 2          4.9         3.0          1.4         0.2 setosa
	# 3          4.7         3.2          1.3         0.2 setosa
	# 4          4.6         3.1          1.5         0.2 setosa
	# 5          5.0         3.6          1.4         0.2 setosa
	# 6          5.4         3.9          1.7         0.4 setosa
	if( !("class"%in%names(dataDf)) ){
		stop("dataDf should have a class column.")
	}
	
	if(is.null(nvars)){ 
		nvars = ncol(dataDf)-1 
	}
	dirDf = data.frame( diag(1,ncol(dataDf)-1,ncol(dataDf)-1) )[,1:nvars]
	names(dirDf) = names(dataDf)[names(dataDf)!="class"][1:nvars]
	plotMatrix2(
			dirDf = dirDf,
			featureDf = dataDf[,names(dataDf)!="class"],
			col = classColors[as.numeric(dataDf$class)],
			alpha=alpha
	)	
	
}
