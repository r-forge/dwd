###############################################################################
## DWD.R
##
## functions to show DWD result
## 
## Author: Xiaosun Lu
###############################################################################

##=============================================================================
## Function to plot data prejection onto the dwd direction and the OPC1
## with two classes. No predictions yet
##=============================================================================

DWDResult = function(
		dataDf, 
		classNames, 
		savePlot=NULL,
		showPlot=TRUE,
#		dwdOption=1, 
		classColors = c("red","blue") ){
	## col = color for two classes
	
	## check the inputs
	if( !("class"%in%names(dataDf)) ){
		stop('dataDf should include class column.')
	}
	if( length(classNames)!=2 ){
		stop('Only consider two classes for now.')
	}
	if( length(classColors)!=length(classNames) ){
		stop('Each class should be assigned a color.')
	}
	if( !unique(classNames%in%dataDf$class) ){
		stop('The class name is not correct.')
	}
	
	
	if(!is.null(savePlot)){ showPlot=TRUE }
	
	
	#############################
	if( length(classNames)==2 ){
		dataDf = subset(dataDf, class%in%classNames)
		dataDf$class = factor(dataDf$class,levels=classNames)
		
		class1Df = subset(dataDf, class==classNames[1])
		class2Df = subset(dataDf, class==classNames[2])
		
		class1Mat = as.matrix(class1Df[,names(class1Df)!="class"])
		class2Mat = as.matrix(class2Df[,names(class2Df)!="class"])
		
#		dim(class1Mat)
#		dim(class2Mat)
		
		#############################################
		## DWD direction
		
#		## option=1: Jason's
#		temp = system.time(
##				output <- dwd(t(class1Mat),t(class2Mat),option=dwdOption)
#				output <- dwd(t(class1Mat),t(class2Mat),option=1)
#		)
		
		temp = system.time(
				output <- kdwd(class~.,data=dataDf)
		)
		
		time = temp[3]
		attr(time,"names")=""
		
		## w is the normal direction vector, beta is the interception term for DWD separation plane
		v = as.vector(output@w[[1]])
		beta <- as.numeric(output@b0)
		dwdDir = -v
		
		#############################################
		## Prediction
		
		dataMat = as.matrix(dataDf[,names(dataDf)!= "class"])
		y = as.vector(dataMat%*%v + as.numeric(beta))
		pred = y
		pred[y>0] = classNames[1]	
		pred[y<0] = classNames[2]
		
		predErrorRate = sum(pred!=dataDf$class)/length(pred)
		
		#############################################
		## Project onto  DWD direction and OPC1 
		
		if(!is.null(savePlot)){
			pdf(savePlot)
		}
		dirDf = NULL
		if(showPlot){
			dirDf = orthoPC(dataMat = dataMat, dir = dwdDir)
			plotMatrix2(
					dirDf = data.frame(OPC1 = dirDf$OPC1, dwdDir = dirDf$dir),
					featureDf = as.data.frame(scale(dataMat,center=TRUE, scale=FALSE)),
					col = classColors[as.numeric(dataDf$class)],
					alpha=1
			)	
		}
		if(!is.null(savePlot)){
			dev.off()
		}
		
		return(
				list(dwdDir = dwdDir, opcDirDf = dirDf[,-1], dwdRunningTime = time, 
						dimClass1 = dim(class1Mat), dimClass2 = dim(class2Mat),
						pred = pred, predErrorRate = predErrorRate,beta=beta)
		)
	}
	
}



##=============================================================================
## DWD predict function
##=============================================================================

## return testDf containing column predClass
DWD.pred = function(trainDf, testDf,dwdOption=0,
		showPlot=FALSE,
		classColors = NULL,testClassShapes=NULL,alphaTrain=NULL){
	## The plot: project all of the data on to dwd direction
	## colored by real classes 
	## The predicted classes of test data is shown by point shape
	## the color of training data is lighter (alphaTrain)
	
#	browser()
	if( !("class"%in%names(trainDf)) ){
		stop('trainDf should include class column.')
	}
	
	classes=levels(trainDf$class)
	if( length(classes)==2 ){
		##############
		## DWD
		class1Df = subset(trainDf, class==classes[1])
		class2Df = subset(trainDf, class==classes[2])
		
		class1Mat = as.matrix(class1Df[,names(class1Df)!="class"])
		class2Mat = as.matrix(class2Df[,names(class2Df)!="class"])
		
		temp = system.time(
				output <- dwd(t(class1Mat),t(class2Mat),option=dwdOption)
		)
		time = temp[3]
		attr(time,"names")=""
		v = as.vector(output$dirvec)
		beta <- as.numeric(output$beta)
		dwdDir = -v
		
		###############
		## predict
		testMat = as.matrix(testDf[,names(testDf)!= "class"])
		y = as.vector(testMat%*%v + as.numeric(beta))
		pred = y
		pred[y>0] = classes[1]	
		pred[y<0] = classes[2]
		
		predErrorRate = sum(pred!=as.character(testDf$class))/length(pred)
		
		## plots
		if(showPlot){
			dDf = rbind(cbind(trainDf,type="train"),cbind(testDf,type="test"))
			dDf[,!names(dDf)%in%c("class","type")] = scale(dDf[,!names(dDf)%in%c("class","type")],center=TRUE,scale=FALSE)
			col1 = classColors[as.numeric(dDf$class)]
			shape1 = rep(1,nrow(dDf))
			shape1[dDf$type=="test"] = testClassShapes[as.numeric(factor(pred, levels=levels(trainDf$class)))]
			alpha1 = rep(1,nrow(dDf))
			alpha1[dDf$type=="train"] = alphaTrain
			dirDf = orthoPC(dataMat = as.matrix(dDf[,!names(dDf)%in%c("class","type")]), dir = dwdDir)
			plotMatrix2(
					dirDf = data.frame(OPC1 = dirDf$OPC1, dwdDir = dirDf$dir),
					featureDf = dDf[,!names(dDf)%in%c("class","type")],
					col = col1,
					shape=shape1,
					alpha=alpha1
			)	
		}#
		
		
	}
	
	return(list(pred=pred,dwdRunningTime = time, dwdDir = dwdDir, beta=beta,
					predErrorRate = predErrorRate))
}

##=============================================================================
## CV to compute misclassification rate using the function above
##=============================================================================
DWD.CVError = function(
		
		dataDf, 
		dwdOption=0, 
		ptrain = 0.8, ## using 80% data to train
		nrep = 10){
	
#	browser()
	
	cvChunkFun = function(i,ddf=dataDf,pt=ptrain,dopt=dwdOption){
		## sampling
		tDf = ddply(ddf, .(class), function(df,pt1=pt){
					size = floor(nrow(df)*pt1)
					data.frame(train= (1:nrow(df)%in%sample(nrow(df), size)) )
				})
		
		trainDf = ddf[tDf$train,]
		testDf = ddf[!tDf$train,]
		
		## dwd
		temp = DWD.pred(trainDf=trainDf, testDf=testDf,dwdOption=dwdOption)
		
		c(
				predErrorRate=temp$predErrorRate,
				dwdRunningTime=temp$dwdRunningTime
		)
	}
	#cvChunkFun(1)
	
	rs = sfSapply(1:nrep,cvChunkFun,ddf=dataDf,pt=ptrain,dopt=dwdOption)
	
	predErrorRate = sum(rs["predErrorRate",])/nrep
	dwdRunningTime = sum(rs["dwdRunningTime",])/nrep
	
	return(list(predErrorRate=predErrorRate,dwdRunningTime=dwdRunningTime,
					predErrorRates = rs["predErrorRate",],dwdRunningTimes = rs["dwdRunningTime",]))
	
}
