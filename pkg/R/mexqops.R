mexqops = function(blk,x,y,options){
  numblk <- length(blk$size)
  cumblk <- c(0,cumsum(blk$size))
  if(options<3)
    nz <- numblk
  else
    nz <- sum(blk$size)
  Z <- .C("mexqops_c",as.double(x),as.double(y),z=double(nz),
          as.integer(numblk),as.integer(cumblk),as.integer(options))
  return(Z$z)
}
