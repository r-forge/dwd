###############################################################################
## simuExample2.R
##
## DWD aplication on simulated two class example
## and compare with other classification tools such as linear SVM and LDA
## 
## Author: Xiaosun Lu
###############################################################################

rm(list=ls())
gc()
library(ggplot2)
library(kernlab)
library(RUnit)
library(DWDnew)
#loadLibraries()
source("sourceLibrary.R")
sessionInfo()
# R version 2.13.0 Patched (2011-06-19 r56186)
# Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
# 
# locale:
# [1] C/C/en_US/C/C/C
# 
# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
#  [1] randomForest_4.6-2    MASS_7.3-13           kernlab_0.9-12       
#  [4] snowfall_1.84         snow_0.3-5            RUnit_0.4.26         
#  [7] Matrix_0.999375-50    lattice_0.19-26       fields_6.3           
# [10] spam_0.23-0           preprocessCore_1.14.0 nortest_1.0          
# [13] outliers_0.14         tseries_0.10-25       zoo_1.6-5            
# [16] quadprog_1.5-4        DWDstable_0.10        Cairo_1.4-9          
# [19] ggplot2_0.8.9         proto_0.3-9.2         reshape_0.8.4        
# [22] plyr_1.5.2            rj_0.5.2-1           
# 
# loaded via a namespace (and not attached):
# [1] rJava_0.8-8  tools_2.13.0
# Warning messages:
# 1: 'DESCRIPTION' file has 'Encoding' field and re-encoding is not possible 
# 2: 'DESCRIPTION' file has 'Encoding' field and re-encoding is not possible 


# CairoX11()
ls()

## I would recommend setting up the parallel computing here
## with any options
## start parallel computing 
## setup the parallel computing outside of the function
nclusters = 15
startCluster(nclusters=nclusters) 

##=============================================================================
## Simulate the data and show plots
##=============================================================================

## nreps = the sample sizes in the two groups
## d = the the number of columns (variables)
## muDist = 1/2* distance between two groups
## classNames = names to be used for the two classes

## Generate a simulated two class data problem
## based on the model described above
classNames = c("N","P")
dataDf = simTwoClassData(nreps=c(100, 100), d=1000, muDist=0.08*sqrt(1000), cNames=classNames,setSeed=33)
dim(dataDf)

0.08*sqrt(1000)
# [1] 2.529822


## the levels of the class variable are the same as the order given in classNames
## the classColors are always applied based on the level of class
classColors = c("red","blue")

## view the data
pdf("output/simu2_originPlotMatrix.pdf")
originPlotMatrix(
		dataDf = dataDf, 
		classColors=classColors,
		nvars=3)
dev.off()


pdf("output/simu2_pcaPlotMatrix.pdf")
pcaPlotMatrix(
		dataDf  = dataDf, 
		classColors=classColors,
		npc=3
)
dev.off()

##==============================================================
## DWD
##==============================================================

## nrepeat = the number of repeats of the CV
nrepeat = 100

## Jason's
rsOpt1.list = DWDResult(dataDf=dataDf,classNames=classNames,
#		dwdOption=1,
		classColors=classColors,
#		showPlot=TRUE)
		savePlot = "output/simu2_1.pdf")

# temp = DWD.CVError(dataDf=dataDf,ptrain=0.8,dwdOption=1,nrep=nrepeat)

###############
## Hanwen's
rsOpt0.list = DWDResult(dataDf=dataDf,classNames=classNames,
#		dwdOption = 0,
		classColors=classColors,
#		showPlot=TRUE)
		savePlot = "output/simu2_0.pdf")
# temp = DWD.CVError(dataDf,dwdOption=0,ptrain=0.8,nrep=nrepeat)

rsOpt1.list$dwdRunningTime
rsOpt0.list$dwdRunningTime
#        
# 43.353 
#       
# 6.195 


## angle between dwd directions of Jason's and Hanwen's
#angle(rsOpt0.list$dwdDir, rsOpt1.list$dwdDir, lessThan90=TRUE)
# [1] 1.815148e-05



##==============================================================
## SVM
##==============================================================

rsSVM.list = lsvmResult(dataDf=dataDf, classNames=classNames, 
		classColors=classColors,cross.svm=0,
#		showPlot=TRUE)
		savePlot = "output/simu2_svm.pdf")
rsSVM.list$svmRunningTime 
#       
# 4.122 

## angle between svm direction and Hanwen's dwd direction
angle(rsSVM.list$svmDir, rsOpt0.list$dwdDir, lessThan90=TRUE)
# [1] 29.93996

######################
### 5-fold CV by ksvm 
#
#temp = lsvmResult(dataDf=dataDf, classNames=classNames, 
#		classColors=classColors,cross.svm=5)
#attr(temp$svmfit,"cross") ## error
## [1] 0.44


##==============================================================
## LDA
##==============================================================

rsLDA.list = ldaResult(dataDf=dataDf,classNames=classNames, 
		classColors=classColors,
#		showPlot=TRUE)
		savePlot = "output/simu2_lda.pdf")
rsLDA.list$ldaRunningTime
#       
# 5.693 


## angle between svm direction and Hanwen's dwd direction
angle(rsLDA.list$ldaDir, rsOpt0.list$dwdDir, lessThan90=TRUE)


rsLDA.list$predErrorRate
# [1] 0.33

## CV
# temp = LDA.CVError(dataDf=dataDf,ptrain=0.8, nrep=nrepeat)
# mean(temp$predErrorRates)
# [1] 0.283

##==============================================================
## compare SVM, DWD and LDA predictions
##==============================================================

## an example of svm prediction to check if my way of computing the 
## prediction (svmPredOpt=1) is the same with predict.ksvm so that 
## the svmDir that I get is correct

## randomly choose training data (80%)
set.seed(10)
## get 80% training data from each class
trainsample = (1:200)%in%c(sample(1:100, 80),sample(101:200,80))
testsample = !trainsample

## view the first dimension
predPlot(trainDf=dataDf[trainsample,],testDf=dataDf[testsample,],
		dir = -c(1,rep(0,ncol(dataDf)-2)),## the first dimension
		pred = dataDf[testsample,]$class,dirName = "dimension_1",
		classColors=classColors, testClassShapes=c(2,3),alphaTrain=0.999,
		savePlot = "output/simu2_predExample.pdf")

## predict.ksvm
temp0 = lsvm.pred(trainDf=dataDf[trainsample,],testDf=dataDf[testsample,],svmPredOpt=0,
		classColors=classColors, testClassShapes=c(2,3),alphaTrain=0.999,
		showPlot=TRUE)
#       showPlot=FALSE)
## my prediction
pdf("output/simu2_svmPred.pdf")
temp1 = lsvm.pred(trainDf=dataDf[trainsample,],testDf=dataDf[testsample,],svmPredOpt=1,
		classColors=classColors, testClassShapes=c(2,3),alphaTrain=0.999,
		showPlot=TRUE)
dev.off()
checkEquals(temp1$predErrorRate,temp0$predErrorRate)

## The prediction error is high.
## One can see the plot above for more details
predsvmError = temp1$predErrorRate



## dwd prediction of the same example above
pdf("output/simu2_dwdPred.pdf")
temp = DWD.pred(trainDf=dataDf[trainsample,],testDf=dataDf[testsample,], dwdOpt=0,
		classColors=classColors, testClassShapes=c(2,3),alphaTrain=0.999,
		showPlot=TRUE)
dev.off()
preddwdError = temp$predErrorRate

## lda prediction of the same example above
pdf("output/simu2_ldaPred.pdf")
temp = lda.pred(trainDf=dataDf[trainsample,],testDf=dataDf[testsample,],
		classColors=classColors, testClassShapes=c(2,3),alphaTrain=0.999,
		showPlot=TRUE)
dev.off()
predldaError = temp$predErrorRate

##=============================================================================
## Compare different classification approaches via CV
##=============================================================================

compare.list = compareFun(dataDf = dataDf,nrep=nrepeat,ptrain=0.8)
names(compare.list)
#  [1] "predErrorRate.j.CI"   "predErrorRate.h.CI"   "predErrorRate.lda.CI"
#  [4] "predErrorRate.svm.CI" "angleDiff.j.CI"       "angleDiff.lda.CI"    
#  [7] "angleDiff.svm.CI"     "predErrorRate.j"      "predErrorRate.h"     
# [10] "predErrorRate.lda"    "predErrorRate.svm"    "angleDiff.j"         
# [13] "angleDiff.lda"        "angleDiff.svm"       


## plot
pdf("output/simu2_cv.pdf")
## plot angle diff between Jason's and Hanwen's dwd
angleDiff = compare.list$angleDiff.j
tempDf = data.frame(x=angleDiff,y=height(angleDiff))
xmax = max(angleDiff)+diff(range(angleDiff))/5
g = ggplot(data=tempDf) + geom_density(aes(x))+ 
		geom_point(aes(x=x,y=y,colour=I("blue"))) +
		geom_vline(xintercept=0,colour="red")+
		geom_vline(xintercept = quantile(angleDiff,0.025),colour="cyan")+
		geom_vline(xintercept = quantile(angleDiff,0.975),colour="cyan")+
		scale_colour_identity()+
		opts(legend.position = "none")+ xlim(c(0,xmax))+
		xlab("Angle between Jason's and Hawen's DWD directions")+ ylab("Density") +
		opts(title=(paste("nrep = ",nrepeat, ", ptrain = ", 0.8,sep="")))
print(g)

order = order(compare.list$predErrorRate.h)
tempDf = data.frame(
		x=rep(1:nrepeat,3),
		y=c(compare.list$angleDiff.j[order],compare.list$predErrorRate.j[order],
				compare.list$predErrorRate.h[order]),
		col=rep(c("red","blue","darkgreen"),each=nrepeat),
		title=rep(c("Angle Diff","Jason's Pred Error Rate","Hanwen's Pred Error Rate"),each=nrepeat)
)
g = ggplot(tempDf) + geom_point(aes(x,y))+ facet_wrap(~title,ncol=1,scales="free")+
		opts(title="Compare Hanwen's DWD results with Jason's")+xlab("CV Replication")+ylab("")
print(g)

dev.off()

##==============================================================
## random forest
##==============================================================
temp=system.time(
		fit.rf <- randomForest(class~.,data=dataDf, ntree=500)
)
rfRunningTime = temp[3]
rfRunningTime

rfPredErrorRate = sum(fit.rf$predicted!=dataDf$class)/nrow(dataDf)
print(fit.rf)

##==============================================================
## save
##==============================================================
## save workspace
save.image(file="rdata/simu2Workspace.RData")


sfStop()
