###############################################################################
## simuData.R
##
## Simulate example data
## 
## Author: Xiaosun Lu
###############################################################################

## I highly recommend that you don't use n as a variable name
## nreps

simTwoClassData = function(nreps,d,
		setSeed = NULL,
		muDist = 1/8*sqrt(d), ## control the distance between two groups,
		cNames=c("P","N")){
	##  nreps is a vector with group sample sizes
	##  d = number of dimensions
	##  muDist = 1/2*group distance
	if( length(nreps)!=2 ){
		stop('Only two-class case is considered for now.')
	}
	
	if(!is.null(setSeed)){set.seed(setSeed)}
	if(length(nreps)==2){
		n1 = nreps[1]
		n2 = nreps[2]
		
		## generate values from a random normal
		datap <- matrix(rnorm(d*n1),ncol=d,nrow=n1)
		datan <- matrix(rnorm(d*n2),ncol=d,nrow=n2)
		## set the means in each group based on the signal size
		datap[,1] <- datap[,1]+rep(muDist,n1)
		datan[,1] <- datan[,1]-rep(muDist,n2)
		
		dataDf <- as.data.frame(rbind(datap,datan));
		## all class column
		dataDf$class = factor(c(rep(cNames[1],n1),rep(cNames[2],n2)),levels=cNames)
		dim(dataDf)
	}
	return(dataDf)
}
#dataDf = simuData(n=c(100, 100), d=1000, muDist=1/8, setSeed=452)









