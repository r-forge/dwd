###############################################################################
## cv.R
##
## Cross validation to compare Jason's DWD, Hanwen's DWD, SVM and LDA
## 
## Author: Xiaosun Lu
###############################################################################

##=============================================================================
## compare differnt classification approaches
##=============================================================================

compareFun = function(
		
		dataDf, 
		ptrain = 0.8, ## using 80% data to train
		nrep = 10
){

    cvCompareChunkFun = function(i,dDf = dataDf, pt=ptrain){
		tDf = ddply(dDf, .(class), function(ddf){
					size = floor(nrow(ddf)*pt)
					data.frame(train= (1:nrow(ddf)%in%sample(nrow(ddf), size)) )
				})
		trainDf = dataDf[tDf$train,]
		testDf = dataDf[!tDf$train,]
		
		## different approaches
		jdwd.list = DWD.pred(trainDf=trainDf, testDf=testDf,dwdOption=1)
		hdwd.list = DWD.pred(trainDf=trainDf, testDf=testDf,dwdOption=0)
		lda.list = lda.pred(trainDf=trainDf, testDf=testDf)
		svm.list = lsvm.pred(trainDf=trainDf, testDf=testDf, svmPredOpt = 1, cross.svm=0,C.svm=1)
		
		jdwdDir = jdwd.list$dwdDir
		hdwdDir = hdwd.list$dwdDir
		ldaDir = lda.list$ldaDir
		svmDir = svm.list$svmDir
		
		
		return(c(
						predErrorRate.j = jdwd.list$predErrorRate,
						predErrorRate.h = hdwd.list$predErrorRate,
						predErrorRate.lda = lda.list$predErrorRate,
						predErrorRate.svm = svm.list$predErrorRate,
						## compare angle diff between each approach and Hanwen's dwd
						angleDiff.j = angle(jdwdDir, hdwdDir, lessThan90=TRUE),
						angleDiff.lda = angle(ldaDir, hdwdDir, lessThan90=TRUE),
						angleDiff.svm = angle(svmDir, hdwdDir, lessThan90=TRUE),
						
						RunningTime.j = jdwd.list$dwdRunningTime,
						RunningTime.h = hdwd.list$dwdRunningTime,
						RunningTime.lda = lda.list$ldaRunningTime,
						RunningTime.svm = svm.list$svmRunningTime
		
				))
	}
	# temp = cvCompareChunkFun(1)
	rs = sfSapply(1:nrep, cvCompareChunkFun,dDf=dataDf,pt=ptrain)

	
	predErrorRate.j = rs["predErrorRate.j",]
	predErrorRate.h = rs["predErrorRate.h",]
	predErrorRate.lda = rs["predErrorRate.lda",]
	predErrorRate.svm = rs["predErrorRate.svm",]
	
	angleDiff.j = rs["angleDiff.j",]
	angleDiff.lda = rs["angleDiff.lda",]
	angleDiff.svm = rs["angleDiff.svm",]
	
	
	return(list(
					predErrorRate.j.CI = c(quantile(predErrorRate.j,0.025),quantile(predErrorRate.j,0.975)),
					predErrorRate.h.CI = c(quantile(predErrorRate.h,0.025),quantile(predErrorRate.h,0.975)),
					predErrorRate.lda.CI = c(quantile(predErrorRate.lda,0.025),quantile(predErrorRate.lda,0.975)),
					predErrorRate.svm.CI = c(quantile(predErrorRate.svm,0.025),quantile(predErrorRate.svm,0.975)),
					
					angleDiff.j.CI = c(quantile(angleDiff.j,0.025),quantile(angleDiff.j,0.975)),
					angleDiff.lda.CI = c(quantile(angleDiff.lda,0.025),quantile(angleDiff.lda,0.975)),
					angleDiff.svm.CI = c(quantile(angleDiff.svm,0.025),quantile(angleDiff.svm,0.975)),
					
					predErrorRate.j = predErrorRate.j,
					predErrorRate.h = predErrorRate.h,
					predErrorRate.lda = predErrorRate.lda,
					predErrorRate.svm = predErrorRate.svm,
					
					angleDiff.j = angleDiff.j,
					angleDiff.lda = angleDiff.lda,
					angleDiff.svm = angleDiff.svm
					
			))
	
}












