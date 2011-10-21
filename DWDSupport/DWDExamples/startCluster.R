##================================================
## startCluster.R
##
## starts the cluster and loads all libraries
## 
## Author: Haaland
##================================================
##Sys.setenv("http_proxy"="http://10.2.96.60")

startCluster = function(nclusters=nclusters){
	sfInit(parallel=TRUE, cpus=nclusters)
	
	sfSource("../DWDtest/setup.R")
	
	sfLibrary(RUnit)
	sfLibrary(plyr)
	## SVM
	sfLibrary(kernlab)
	## LDA
	sfLibrary(MASS)
	
	sfSource("../DWDUtility/sourceLibrary.R")
	
	#sfExportAll()
	
}

loadLibraries = function(){
	
	library(ggplot2)
	library(Cairo)
	
	source("../DWDtest/setup.R")
	
	library(RUnit)
	library(snowfall) ## install.packages("snowfall", dep = TRUE, repos = "http://cran.stat.ucla.edu/")
	
	## SVM
	library(kernlab)
	## LDA
	library(MASS)
	## RF
	library(randomForest)
	
	source("../DWDUtility/sourceLibrary.R")
	
}
