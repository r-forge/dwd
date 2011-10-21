##********************************************************************
## checkdense : identify the dense columns of a matrix
## 
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##********************************************************************

checkdense=function(A){
	m <- dim(A)[1]
	idxden <- NULL
	nzratio <- 1
	if (m > 1000)
		nzratio <- 0.20
	if (m > 2000)
		nzratio <- 0.10
	if (m > 5000)
		nzratio <- 0.05
	if (nzratio < 1){
		nzcolA <- apply(A,2,nnz)
		idxden <- which(nzcolA > nzratio*m)
		if(length(idxden)==0)
			idxden <- NULL
	}
	return(idxden)
}

nnz <- function(x){
	return(sum(abs(x)>0))
}
