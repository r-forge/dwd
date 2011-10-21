###############################################################################
## SVM.R
##
## functions to show SVM results
## 
## Author: Xiaosun Lu
###############################################################################

##=============================================================================
## Function to plot data prejection onto the svm direction and the OPC1
## with two classes. No predictions yet
##=============================================================================


lsvmResult = function(
		dataDf, classNames, classColors = c("red","blue"),
		savePlot = NULL, showPlot = TRUE,
		C.svm=1, cross.svm=0,xlims=NULL){
	
	## xlims = a list of min, and max of x axis in the plot
	
	if( !("class"%in%names(dataDf)) ){
		stop('trainDf should include class column.')
	}
	if( length(classNames)!=2 ){
		stop('Only two-class classification is considered for now.')
	}
	if( length(classColors)!=length(classNames) ){
		stop('Each class should be assigned a color.')
	}
	if( !unique(classNames%in%dataDf$class) ){
		stop('The class name is not correct.')
	}
	if(!is.null(savePlot)){ showPlot=TRUE }
	
	#####################
	
	if( length(classNames)==2 ){
		dDf = subset(dataDf,class%in%classNames)
		dDf$class = factor(dDf$class,levels=classNames)
		# table(dDf$class)
		# 
		# versicolor  virginica 
		#         50         50 
		
		# head(dDf)
		#    Sepal.Length Sepal.Width Petal.Length Petal.Width      class
		# 51          7.0         3.2          4.7         1.4 versicolor
		# 52          6.4         3.2          4.5         1.5 versicolor
		# 53          6.9         3.1          4.9         1.5 versicolor
		# 54          5.5         2.3          4.0         1.3 versicolor
		# 55          6.5         2.8          4.6         1.5 versicolor
		# 56          5.7         2.8          4.5         1.3 versicolor
		
		##=============================================================================
		## Linear SVM C classification on dDf
		##=============================================================================
		
		## A high cost value C will force the SVM to create a complex enough prediction 
		## function to missclassify as few training points as possible, while a lower 
		## cost parameter will lead to a simpler prediction function
		
		## cross = 5 : 5-fold cross validation
		
		##set.seed(setSeed)
		temp = system.time(
				svmfit <- ksvm(class~.,data=dDf,type="C-svc",
						scaled=FALSE,
						kernel="vanilladot",C=C.svm,cross=cross.svm)
		)
		# svmfit
		# Support Vector Machine object of class "ksvm" 
		# 
		# SV type: C-svc  (classification) 
		#  parameter : cost C = 1 
		# 
		# Linear (vanilla) kernel function. 
		# 
		# Number of Support Vectors : 16 
		# 
		# Objective Function Value : -11.2782 
		# Training error : 0.04 
		# Cross validation error : 0.05 
		time = temp[3]
		attr(time,"names")=""
		
		
		##=============================================================================
		## find the lsvm direction
		##=============================================================================
		
		## support vectors
		svDf = (dDf[,names(dDf)!="class"])[attr(svmfit,"alphaindex")[[1]],]
		v = t(attr(svmfit,"coef")[[1]]%*%as.matrix(svDf))
		svmDir = as.vector(v/sqrt(sum(v^2)))
		
		#############################################
		## Project onto  SVM direction and OPC1 
		
		if(!is.null(savePlot)){
			## results
			pdf(savePlot)
		}
		dirDf = NULL
		if(showPlot){
			dirDf = orthoPC(dataMat = as.matrix(dDf[,names(dDf)!="class"]), dir = svmDir)
			plotMatrix2(
					dirDf = data.frame(OPC1 = dirDf$OPC1, svmDir = dirDf$dir),
					featureDf = scale(dDf[,names(dDf)!="class"],center=TRUE, scale=FALSE),
					col = classColors[as.numeric(dDf$class)],
					alpha=1,
					xlims=xlims
			)	
		}
		if(!is.null(savePlot)){
			dev.off()
		}
		
		#########################################
		## prediction
		
		pred = predict(svmfit)
		predErrorRate = sum(pred!=dDf$class)/length(pred)
		
		
		return(list(svmDir = svmDir, opcDirDf = dirDf[,-1], svmRunningTime=time,
						svmfit = svmfit,  pred = pred, predErrorRate = predErrorRate
				))
	}
}

##==============================================================
## SVM prediction
##==============================================================

## opt=0, using predict.ksvm, which seems not right for high dim example 
lsvm.pred = function(trainDf, testDf, svmPredOpt = 0, cross.svm=0,C.svm=1,
		showPlot=FALSE,
		classColors = NULL,testClassShapes=NULL,alphaTrain=NULL){
	## The plot: project all of the data on to svm direction
	## colored by real classes 
	## The predicted classes of test data is shown by point shape
	## the color of training data is lighter (alphaTrain)
	
#	 browser()
	classes=levels(trainDf$class)
	if( length(classes)==2 ){
		temp = system.time(
				fit.svm <- ksvm(class~.,data=trainDf,type="C-svc",kernel="vanilladot",
						scaled = FALSE,
						C=C.svm,cross=cross.svm)
		)
		time = temp[3]
		attr(time,"names")=""
		
		## svmDir
		svDf = (trainDf[,names(trainDf)!="class"])[attr(fit.svm,"alphaindex")[[1]],]
		v = t(attr(fit.svm,"coef")[[1]]%*%as.matrix(svDf))
		svmDir = as.vector(v/sqrt(sum(v^2)))
		
		## predict
		if(svmPredOpt==0){
			pred = predict(fit.svm,testDf)
		}
		if(svmPredOpt==1){
			## predict
			testMat = as.matrix(testDf[,names(testDf)!= "class"])
			y = -as.vector(testMat%*%v - as.numeric(attr(fit.svm,"b")))
			pred = y
			pred[y>0] = classes[1]	
			pred[y<0] = classes[2]
		}
		predErrorRate = sum(pred!=testDf$class)/length(pred)
		
		
		## plots
		if(showPlot){ 
			dDf = rbind(cbind(trainDf,type="train"),cbind(testDf,type="test"))
			dDf[,!names(dDf)%in%c("class","type")] = scale(dDf[,!names(dDf)%in%c("class","type")],center=TRUE,scale=FALSE)
			col1 = classColors[as.numeric(dDf$class)]
			shape1 = rep(1,nrow(dDf))
			shape1[dDf$type=="test"] = testClassShapes[as.numeric(factor(pred, levels=levels(trainDf$class)))]
			alpha1 = rep(1,nrow(dDf))
			alpha1[dDf$type=="train"] = alphaTrain
			dirDf = orthoPC(dataMat = as.matrix(dDf[,!names(dDf)%in%c("class","type")]), dir = svmDir)
			plotMatrix2(
					dirDf = data.frame(OPC1 = dirDf$OPC1, svmDir = dirDf$dir),
					featureDf = dDf[,!names(dDf)%in%c("class","type")],
					col = col1,
					shape=shape1,
					alpha=alpha1
			)	
		}#
		
		
		return(list(pred=pred,svmRunningTime = time, svmDir = svmDir,
						predErrorRate = predErrorRate))
	}
}

##=============================================================================
## LSVM CV error
##=============================================================================
lsvm.CVError = function(
		
		dataDf, 
		ptrain = 0.8, ## using 80% data to train
		nrep = 10,
		svmPredOpt = 0){
	
	cvChunkFun = function(i, dDf = dataDf,pt=ptrain,spOpt = svmPredOpt){
		tDf = ddply(dDf, .(class), function(ddf){
					size = floor(nrow(ddf)*pt)
					data.frame(train= (1:nrow(ddf)%in%sample(nrow(ddf), size)) )
				})
		trainDf = dataDf[tDf$train,]
		testDf = dataDf[!tDf$train,] ## test
		
		temp = lsvm.pred(trainDf = trainDf, testDf = testDf, svmPredOpt = svmPredOpt, cross.svm=0, C.svm=1)
		return(predErrorRate=temp$predErrorRate)
	}
	# cvChunkFun(1)
	errors = sfSapply(1:nrep, cvChunkFun, dDf = dataDf, pt=ptrain,spOpt = svmPredOpt)
	return(list(predErrorRates=errors))
}



