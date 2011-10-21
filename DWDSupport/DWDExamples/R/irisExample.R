###############################################################################
## irisExample.R
##
## DWD application on iris data
## 
## Author: luxi
###############################################################################



rm(list=ls())
gc()

source("../DWDUtility/startCluster.R")
loadLibraries()

## start parallel computing
nclusters = 15
startCluster(nclusters=nclusters) 

#=============================================================================
# Data
#=============================================================================
data(iris)

dataDf = data.frame(iris)
names(dataDf)[5] = "class"
head(dataDf)
#   Sepal.Length Sepal.Width Petal.Length Petal.Width  class
# 1          5.1         3.5          1.4         0.2 setosa
# 2          4.9         3.0          1.4         0.2 setosa
# 3          4.7         3.2          1.3         0.2 setosa
# 4          4.6         3.1          1.5         0.2 setosa
# 5          5.0         3.6          1.4         0.2 setosa
# 6          5.4         3.9          1.7         0.4 setosa

table(dataDf$class)
# 
#     setosa versicolor  virginica 
#         50         50         50 
dataDf$class = factor(dataDf$class,levels=c("setosa", "versicolor", "virginica"))


##==============================================================
## view the data
##==============================================================

pdf("output/iris_originPlotMatrix.pdf")
originPlotMatrix(
		dataDf = dataDf, 
		classColors=c("purple","red","blue")
)
dev.off()


pdf("output/iris_pcaPlotMatrix.pdf")
pcaPlotMatrix(
		dataDf  = dataDf, 
		classColors =c("purple", "red","blue")
)
dev.off()


## Now I only focus on the two classes that seem hardest distinguish
classNames = c("versicolor","virginica")
dataDf = subset(dataDf,class%in%classNames)
dataDf$class = factor(dataDf$class,levels=classNames)

classColors = c("red","blue")




##==============================================================
## DWD
##==============================================================

## Hanwen's
rsOpt0.list = DWDResult(dataDf=dataDf,classNames=classNames,
		classColors=classColors,dwdOption = 0,
		savePlot = "output/iris_VersicolorVirginica_0.pdf")
## Jason's
rsOpt1.list = DWDResult(dataDf=dataDf,classNames=classNames,
		classColors=classColors,dwdOption = 1,
		savePlot = "output/iris_VersicolorVirginica_1.pdf")
## angle between dwd directions of Jason's and Hanwen's
angle(rsOpt0.list$dwdDir, rsOpt1.list$dwdDir, lessThan90=TRUE)


##=============================================================================
## SVM
##=============================================================================

rsSVM.list = lsvmResult(dataDf=dataDf, classNames=classNames,cross.svm=0, 
		classColors = classColors,
		savePlot = "output/iris_VersicolorVirginica_svm.pdf")
rsSVM.list$svmRunningTime 

## angle between svm direction and Hanwen's dwd direction
angle(rsSVM.list$svmDir, rsOpt0.list$dwdDir, lessThan90=TRUE)



##==============================================================
## LDA
##==============================================================


rsLDA.list = ldaResult(dataDf=dataDf, classNames=classNames, classColors = classColors,
		savePlot = "output/iris_VersicolorVirginica_lda.pdf")
rsLDA.list$ldaRunningTime

## angle between svm direction and Hanwen's dwd direction
angle(rsLDA.list$ldaDir, rsOpt0.list$dwdDir, lessThan90=TRUE)

##=============================================================================
## Compare different classification approaches via CV
##=============================================================================

## nrepeat = the number of repeats of the CV
nrepeat = 100

compare.list = compareFun(dataDf = dataDf,nrep=nrepeat,ptrain=0.8)
names(compare.list)
#  [1] "predErrorRate.j.CI"   "predErrorRate.h.CI"   "predErrorRate.lda.CI"
#  [4] "predErrorRate.svm.CI" "angleDiff.j.CI"       "angleDiff.lda.CI"    
#  [7] "angleDiff.svm.CI"     "predErrorRate.j"      "predErrorRate.h"     
# [10] "predErrorRate.lda"    "predErrorRate.svm"    "angleDiff.j"         
# [13] "angleDiff.lda"        "angleDiff.svm"       


## plot
pdf("output/iris_VersicolorVirginica_cv.pdf")
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
# elapsed 
#   0.029 

rfPredErrorRate = sum(fit.rf$predicted!=dataDf$class)/nrow(dataDf)
print(fit.rf)
# 
# Call:
#  randomForest(formula = class ~ ., data = dataDf, ntree = 500) 
#                Type of random forest: classification
#                      Number of trees: 500
# No. of variables tried at each split: 2
# 
#         OOB estimate of  error rate: 6%
# Confusion matrix:
#            versicolor virginica class.error
# versicolor         47         3        0.06
# virginica           3        47        0.06

##==============================================================
## compare SVM, DWD and LDA predictions
##==============================================================

## an example of svm prediction to check if my way of computing the 
## prediction (svmPredOpt=1) is the same with predict.ksvm so that 
## the svmDir that I get is correct

## randomly choose training data (80%)
set.seed(10)
## get 80% training data from each class
trainsample = (1:100)%in%c(sample(1:50, 40),sample(51:100,40))
testsample = !trainsample
## predict.ksvm
temp0 = lsvm.pred(trainDf=dataDf[trainsample,],testDf=dataDf[testsample,],svmPredOpt=0,
		classColors=classColors, testClassShapes=c(2,3),alphaTrain=0.999,
		showPlot=TRUE)
#       showPlot=FALSE)
## my prediction
pdf("output/iris_svmPred.pdf")
temp1 = lsvm.pred(trainDf=dataDf[trainsample,],testDf=dataDf[testsample,],svmPredOpt=1,
		classColors=classColors, testClassShapes=c(2,3),alphaTrain=0.999,
		showPlot=TRUE)
dev.off()
checkEquals(temp1$predErrorRate,temp0$predErrorRate)

## The prediction error is high.
## One can see the plot above for more details
predsvmError = temp1$predErrorRate



## dwd prediction of the same example above
pdf("output/iris_dwdPred.pdf")
temp = DWD.pred(trainDf=dataDf[trainsample,],testDf=dataDf[testsample,], dwdOpt=0,
		classColors=classColors, testClassShapes=c(2,3),alphaTrain=0.999,
		showPlot=TRUE)
dev.off()
preddwdError = temp$predErrorRate

## lda prediction of the same example above
pdf("output/iris_ldaPred.pdf")
temp = lda.pred(trainDf=dataDf[trainsample,],testDf=dataDf[testsample,],
		classColors=classColors, testClassShapes=c(2,3),alphaTrain=0.999,
		showPlot=TRUE)
dev.off()
predldaError = temp$predErrorRate


##==============================================================
## save
##==============================================================
## save workspace

save.image(file="rdata/irisWorkspace.RData")

sfStop()




