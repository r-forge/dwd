##***********************************************************************
## validate: validate data
##
## [At,C,dim,numblk,X0,Z0] = validate(blk,At,C,b,X0,y0,Z0);
## 
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##***********************************************************************

validate=function(blk,At,C,b,X0,y0,Z0,spdensity){
	
	if(is.null(spdensity))
		spdensity <- 0.5
	
	m <- length(b)  
	for(p in 1:length(blk$type)){
		n <- sum(blk$size[[p]]) 
		if((blk$type[p]=="q") | (blk$type[p]=="l") | (blk$type[p]=="u")){
			if (length(C[[p]])[1] != n) 
				print("validate: blk and C are not compatible") 
			if(dim(At[[p]])[1]==m & dim(At[[p]])[2]==n & m!=n) 
				At[[p]] <- t(At[[p]])
			if(!all(dim(At[[p]]) == c(n,m))) 
				print("validate: blk and At not compatible")
			if(!is(At[[p]],"sparseMatrix"))
				At[[p]] <- as(At[[p]],"dgTMatrix")
#			print(C[[p]])
			if (nnzero(C[[p]]) < spdensity*n){ 
				if(!is(C[[p]],"sparseMatrix"))
					C[[p]] <- as(C[[p]],"sparseMatrix")
			}else{
				if(is(C[[p]],"sparseMatrix"))
					C[[p]] <- as(C[[p]],"denseMatrix")
			}
			if (nargs()==8){
				if(!all(c(dim(X0[[p]])[2],dim(Z0[[p]])[2])==1))
					print(paste("validate: ",as.character(p),
									"-th block of X0,Z0 must be column vectors",sep=""))
				if(!all(c(dim(X0[[p]])[1],dim(Z0[[p]])[1])==n))
					print("validate: blk, and X0,Z0, are not compatible")
				if (nnzero(X0[[p]]) < spdensity*n){ 
					if(class(X0[[p]])!="gdCMatrix")
						X0[[p]] <- as(X0[[p]],"sparseMatrix")
				}else{
					if(class(X0[[p]])=="gdCMatrix")
						X0[[p]] <- as(X0[[p]],"denseMatrix")
				}
				if (nnzero(Z0[[p]]) < spdensity*n){ 
					if(class(X0[[p]])!="gdCMatrix")
						Z0[[p]] <- as(Z0[[p]],"sparseMatrix")
				}else{
					if(class(Z0[[p]])=="gdCMatrix")
						Z0[[p]] <- as(Z0[[p]],"denseMatrix")
				}
				if(blk$type[p]=="u") 
					Z0[[p]] <- as(zeros(n,1),"sparseMatrix")
			}
		}else{
			print(" blk: some fields are not specified correctly")
		}
	}
	
	##-----------------------------------------
	## problem dimension
	##-----------------------------------------
	
	nn <- rep(0,length(blk$type))
	dim <- zeros(1,3)
	numblk <- 0  
	for(p in 1:length(blk$type)){
		if(blk$type[p]=="q"){
			dim[1] <- dim[1] + sum(blk$size[[p]]) 
			numblk <- numblk + length(blk$size[[p]])
			nn[p] <- length(blk$size[[p]])
		}else if(blk$type[p]=="l"){
			dim[2] <- dim[2] + sum(blk$size[[p]]) 
			nn[p] <- sum(blk$size[[p]]) 
		}else if(blk$type[p]=="u"){
			dim[3] <- dim[3] + sum(blk$size[[p]]) 
			nn[p] <- sum(blk$size[[p]]) 
		}
	}
	if(nargs()==5)
		return(list(At=At,C=C,dim=dim,numblk=numblk,
                            spdensity=spdensity))
	else if(nargs()==8)
		return(list(At=At,C=C,dim=dim,numblk=numblk,
                            X=X0,Z=Z0,spdensity=spdensity))
}
