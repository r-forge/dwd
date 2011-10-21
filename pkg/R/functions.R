ones = function(m,n){
	return(matrix(1,nrow=m,ncol=n))
}

zeros = function(m,n){
	return(matrix(0,nrow=m,ncol=n))
}

speye = function(m){
#This function is not actually sparse yet.  Needs to be fixed.
	
	return(diag(m))
	
}
