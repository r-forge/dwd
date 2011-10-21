###############################################################################
## LDA.R
##
## function to show linear discriminant analysis results
## 
## Author: Xiaosun Lu
###############################################################################





ldaResult = function(
		dataDf,classNames,classColors = c("red","blue"),
		savePlot = NULL, showPlot = TRUE,
		option = 1

){
	
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
	if(!is.null(savePlot)){ showPlot=TRUE }
	
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
		
		
		
		## use MASS package to get ldaDir
		temp = system.time(
				fit <- lda(class~.,dDf)
		)
		d = fit$scaling
		d = d/sqrt(sum(d^2))
		ldaDir = as.vector(d)
		
		time = temp[3]
		attr(time,"names")=""
		
#	## another way to compute ldaDir, same result as above when nsample < dimension
#    r=(1.0e-4)
#	
#	x1  = subset(dDf,class==classes[1])[,names(dDf)!="class"]
#	x2  = subset(dDf,class==classes[2])[,names(dDf)!="class"]
#	n1 = nrow(x1)
#	n2 = nrow(x2)
#	
#	x1ba = matrix(rep(mean(x1),each=n1),n1,ncol(x1))
#	x2ba = matrix(rep(mean(x2),each=n2),n2,ncol(x2))
#	
#	s1 = as.matrix(t(x1-x1ba))%*%as.matrix((x1-x1ba))
#	s2 = as.matrix(t(x2-x2ba))%*%as.matrix((x2-x2ba))
#	
#	s = s1 + s2
#	rmat = diag(rep(r,nrow(s)))
#	l = solve(s+rmat)%*%(mean(x2)-mean(x1))
#	l = l/sqrt(sum(l^2))
#	colnames(l) = "LD1"
#	ldaDir = as.vector(l)
#	
#	angle(l,d,lessThan90=TRUE)
#	checkEquals(l,d)
#	# [1] TRUE
		
		
		#############################################
		## Project onto  LDA direction and OPC1 
		
		
		if(!is.null(savePlot)){
			## results
			pdf(savePlot)
		}
		dirDf = NULL
		if(showPlot){
			dirDf = orthoPC(dataMat = as.matrix(dDf[,names(dDf)!="class"]), dir = ldaDir)
			plotMatrix2(
					dirDf = data.frame(OPC1 = dirDf$OPC1, ldaDir = dirDf$dir),
					featureDf = scale(dDf[,names(dDf)!="class"],center=TRUE, scale=FALSE),
					col = classColors[as.numeric(dDf$class)],
					alpha=1
			)	
		}
		if(!is.null(savePlot)){
			dev.off()
		}
		
		#########################################
		## prediction
		pred = predict(fit)$class
		predErrorRate = sum(pred!=dDf$class)/length(pred)
		
		return(list(ldaDir = ldaDir, opcDirDf = dirDf[,-1], ldaRunningTime=time,
						ldafit = fit, pred = pred, predErrorRate = predErrorRate
				))
	}
}
##==============================================================
## LDA prediction
##==============================================================

 
lda.pred = function(trainDf, testDf,
		showPlot=FALSE,
		classColors = NULL,testClassShapes=NULL,alphaTrain=NULL){
	
#	browser()
	
	classes = as.character(unique(trainDf$class))
	
	## lda
	temp = system.time(
			fit.lda <- lda(class~.,trainDf)
	)
	time = temp[3]
	attr(time,"names")=""
	
	## dir
	d = fit.lda$scaling
	d = d/sqrt(sum(d^2))
	ldaDir = as.vector(d)
	
	## pred
	pred = predict(fit.lda,testDf)$class
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
		dirDf = orthoPC(dataMat = as.matrix(dDf[,!names(dDf)%in%c("class","type")]), dir = ldaDir)
		plotMatrix2(
				dirDf = data.frame(OPC1 = dirDf$OPC1, ldaDir = dirDf$dir),
				featureDf = dDf[,!names(dDf)%in%c("class","type")],
				col = col1,
				shape=shape1,
				alpha=alpha1
		)	
	}#
	
	return(list(pred=pred,ldaRunningTime = time, ldaDir = ldaDir,
					predErrorRate = predErrorRate))
}


##=============================================================================
## CV to compute misclassification rate
##=============================================================================
LDA.CVError = function(
		
		dataDf, 
		ptrain = 0.8, ## using 80% data to train
		nrep = 100){
	
	dataDf$class = factor(dataDf$class)
	
	cvChunkFun = function(i, dDf = dataDf,pt=ptrain){
		tDf = ddply(dDf, .(class), function(ddf){
					size = floor(nrow(ddf)*pt)
					data.frame(train= (1:nrow(ddf)%in%sample(nrow(ddf), size)) )
				})
		trainDf = dataDf[tDf$train,]
		testDf = dataDf[!tDf$train,]
		lda.pred(trainDf = trainDf,testDf = testDf)$predErrorRate
	}
	cvChunkFun(1)
	errors = sfSapply(1:nrep, cvChunkFun, dDf = dataDf, pt=ptrain)
	
	
	return(list(predErrorRates=errors))
	
	
}









