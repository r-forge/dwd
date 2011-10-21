###############################################################################
## simuExample.R
##
## DWD aplication (Hanwen's) on simulated two class example
## in order to show that SVM has data piling problem 
## 
## 
## Author: Xiaosun Lu
###############################################################################
rm(list=ls())
gc()

library(animation)

source("../DWDUtility/startCluster.R")
loadLibraries()

##==============================================================
## function to make a movie showing SVM classification by changing d = dimension 
##==============================================================

##==============================================================
## one needs to install ImageMagick in order to use its
## convert function to create the gif movie. However I have 
## no permission to change the system path (it should probably be
## opt/local/bin in $ sudo vi /etc/paths  ), so I copied "convert" under 
## usr/bin, and it works well for now.
##
## to view the movie, it is good to open the file in firefox
##==============================================================



movie = function(opt){
	
	classNames = c("N","P")
	classColors = c("red","blue")
	nsample = 20
	## 1/2 distance between two groups
	muDist = 2.2
	## bandwidth of density
	bw = 0.5
	d = as.integer(c(10,20,30,40,50,100,200, 300, 400, 600, 800,1000,2000))
	
	#############
	## create data
	
	npic = length(d)
	projList = list()
	for (i in 1:npic) {
		# i=1
		## Simulated Data
		dataDf = simTwoClassData(nreps=c(nsample,nsample),setSeed=NULL,
				muDist=muDist,d=d[i],cNames=classNames)
		if(opt=="SVM"){
			dir = lsvmResult(dataDf=dataDf, cross.svm=0,C.svm=1,showPlot=FALSE,
					classNames=classNames)$svmDir
		}else if(opt=="DWD"){
			dir = DWDResult(dataDf=dataDf, dwdOption=0,showPlot=FALSE,
					classNames=classNames)$dwdDir	
		}else if(opt=="LDA"){
			dir = ldaResult(dataDf=dataDf,showPlot=FALSE,classNames=classNames)$ldaDir
		}
		## proj
		dataMat = scale(as.matrix(dataDf[,names(dataDf)!="class"]),center=TRUE, scale=FALSE)
		proj = as.vector(dataMat%*%dir)
		projList[[i]] = proj		
	}
	limMat = sapply(projList,function(proj){
				c(xmin = min(proj)-mean(proj),
						xmax = max(proj)-mean(proj), 
						ymax = max(density(proj,bw=bw)$y))
			})
	xmin = min(limMat[1,])-0.001
	xmax = max(limMat[2,])+0.001
	ymax = max(limMat[3,])+0.001
	
	##############
	## plot
	for (i in 1:npic) {
		proj = projList[[i]]
		temp1Df = data.frame(x=proj,y=runif(nsample*2,ymax/4,ymax/4*2))
		col =  classColors[c(rep(1,nsample),rep(2,nsample))]
		temp1Df = cbind(temp1Df,col)
		dx = density(proj,bw=bw)$x
		dy = density(proj,bw=bw)$y
		temp2Df = data.frame(dx,dy)[dx>xmin & dx<xmax,]
		g = ggplot(temp1Df,aes(x=x,y=y)) + geom_line(data=temp2Df,aes(dx,dy),alpha=0.7)+ 
				geom_point(aes(x=x,y=y,colour=I(col))) +
				scale_colour_identity()+xlim(c(xmin,xmax))+
				opts(legend.position = "none")+ xlab(paste(opt,"Dir",sep=""))+ ylab("density")+
				ylim(c(0,ymax)) + opts(title=paste(opt,"Classification: dim =", d[i]))
		print(g)
	}
	invisible(NULL)
}

saveMovie(movie(opt="SVM"),interval = 0.6, 
		outdir = paste(getwd(),"/output",sep=""),
		movie.name="simu_movie_svm.gif"
)
saveMovie(movie(opt="DWD"),interval = 0.6, 
		outdir = paste(getwd(),"/output",sep=""),
		movie.name="simu_movie_dwd.gif"
)
saveMovie(movie(opt="LDA"),interval = 0.6, 
		outdir = paste(getwd(),"/output",sep=""),
		movie.name="simu_movie_lda.gif"
)

##==============================================================
## plot all three classifications above in one movie
##==============================================================


movieCompare = function(){
	
	classNames = c("N","P")
	classColors = c("red","blue")
	## nsample = number of samples from each group
	nsample = 20
	## muDist = 1/2 distance between two groups
	muDist = 2.2
	## bandwidth of density
	bw = 0.5
	d = as.integer(c(10,20,30,40,50,100,200, 300, 400, 600, 800,1000,2000))
	
	#############
	## create data
	
	npic = length(d)
	dprojList = list()
	sprojList = list()
	lprojList = list()
	for (i in 1:npic) {
		# i=1
		## Simulated Data
		dataDf = simTwoClassData(nreps=c(nsample,nsample),setSeed=NULL,
				muDist=muDist,d=d[i],cNames=classNames)
		## svm dwd and lda directions
		sdir = lsvmResult(dataDf=dataDf, cross.svm=0,C.svm=1,showPlot=FALSE,
				classNames=classNames)$svmDir
		ddir = DWDResult(dataDf=dataDf, dwdOption=0,showPlot=FALSE,
				classNames=classNames)$dwdDir	
		ldir = ldaResult(dataDf=dataDf,showPlot=FALSE,classNames=classNames)$ldaDir
		## proj
		dataMat = scale(as.matrix(dataDf[,names(dataDf)!="class"]),center=TRUE, scale=FALSE)
		sprojList[[i]] = as.vector(dataMat%*%sdir)	
		dprojList[[i]] = as.vector(dataMat%*%ddir)	
		lprojList[[i]] = as.vector(dataMat%*%ldir)	
		
	}
	limMat = sapply(c(sprojList,dprojList,lprojList),function(proj){
				c(xmin = min(proj)-mean(proj),
						xmax = max(proj)-mean(proj), 
						ymax = max(density(proj,bw=bw)$y))
			})
	xmin = min(limMat[1,])+0.001
	xmax = max(limMat[2,])+0.001
	ymax = max(limMat[3,])+0.001
	
	##############
	## plot
	for (i in 1:npic) {
		sproj = sprojList[[i]]
		dproj = dprojList[[i]]
		lproj = lprojList[[i]]
		
		temp1Df = data.frame(x=c(sproj,dproj,lproj),y=rep(runif(nsample*2,ymax/4,ymax/4*2),3),
				type=rep(c("SVM","DWD","LDA"),each=nsample*2))
		temp1Df$type = factor(temp1Df$type,levels=c("DWD","SVM","LDA"))
		col =  rep(classColors[c(rep(1,nsample),rep(2,nsample))],3)
		temp1Df = cbind(temp1Df,col)
		
		
		dx = c(density(sproj,bw=bw)$x,density(dproj,bw=bw)$x,density(lproj,bw=bw)$x)
		dy = c(density(sproj,bw=bw)$y,density(dproj,bw=bw)$y,density(lproj,bw=bw)$y)
		type = rep(c("SVM","DWD","LDA"), each=512)
		temp2Df = data.frame(dx,dy,type)[dx>xmin & dx<xmax,]
		temp2Df$type = factor(temp2Df$type,levels=c("DWD","SVM","LDA"))
		
		
		g = ggplot(temp1Df) + geom_line(data=temp2Df,aes(x=dx, y=dy),alpha=0.7)+ 
				geom_point(aes(x=x,y=y,colour=I(col))) +
				facet_wrap(~type,ncol=1)+
				scale_colour_identity()+xlim(c(xmin,xmax))+
				opts(legend.position = "none")+ xlab("separating directions")+ ylab("density")+
				ylim(c(0,ymax)) + opts(title=paste("Compare Three Classifications: dim =", d[i]))
		print(g)
	}
	invisible(NULL)
}

saveMovie(movieCompare(),interval = 0.6, 
		outdir = paste(getwd(),"/output",sep=""),
		movie.name="simu_movie.gif",
		ani.width=1000, ani.height=800
)

















