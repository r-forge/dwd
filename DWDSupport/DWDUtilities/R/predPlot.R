###############################################################################
## predPlot.R
##
## function to plot training data and test data projection on 
## the separating direction and OPC1
##
## It is the same code that is used in svm.pred, dwd.pred and lda.pred
## 
# Author: luxi
###############################################################################



predPlot = function(trainDf, testDf, dir,savePlot=NULL,pred,
		classColors, testClassShapes,alphaTrain=0.999,dirName){
	
#	## example
#	trainDf=dataDf[trainsample,]
#	testDf=dataDf[testsample,]
#	dir = -c(1,rep(0,ncol(dataDf)-2))
#	pred = dataDf[testsample,]$class
#	dirName = "dimension_1"
#	classColors=classColors
#	testClassShapes=c(2,3)
#	alphaTrain=0.999
#	savePlot = "output/simu2_predExample.pdf"
	
	
	
	
#	browser()
	checkEquals(length(pred),nrow(testDf))
	
	if(!is.null(savePlot)){
		pdf(savePlot)
	}
	dDf = rbind(cbind(trainDf,type="train"),cbind(testDf,type="test"))
	dDf[,!names(dDf)%in%c("class","type")] = scale(dDf[,!names(dDf)%in%c("class","type")],center=TRUE,scale=FALSE)
	col1 = classColors[as.numeric(dDf$class)]
	shape1 = rep(1,nrow(dDf))
	shape1[dDf$type=="test"] = testClassShapes[as.numeric(factor(pred, levels=levels(trainDf$class)))]
	alpha1 = rep(1,nrow(dDf))
	alpha1[dDf$type=="train"] = alphaTrain
	dirDf = orthoPC(dataMat = as.matrix(dDf[,!names(dDf)%in%c("class","type")]), dir = dir)
	dirDf = data.frame(OPC1 = dirDf$OPC1, dir = dirDf$dir)
	names(dirDf)[2] = dirName
	plotMatrix2(
			dirDf = dirDf,
			featureDf = dDf[,!names(dDf)%in%c("class","type")],
			col = col1,
			shape=shape1,
			alpha=alpha1
	)	
	if(!is.null(savePlot)){
		dev.off()
	}
	
}

