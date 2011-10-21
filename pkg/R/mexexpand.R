mexexpand = function(blksize,x){
  numblk <- length(blksize)
  Z <- .C("mexexpand_c",as.integer(blksize),as.integer(numblk),
          as.double(x),z=double(sum(blksize)))
  return(Z$z)
}
