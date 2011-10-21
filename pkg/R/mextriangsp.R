mextriangsp = function(U,b,options){
  n <- length(b)
  U <- as(U,"dgCMatrix")
  Z <- .C("mextriangsp_c",as.double(U@x),as.integer(U@i),as.integer(U@p),
          as.double(b),x=double(n),as.integer(n),
          as.integer(options))
  return(Z$x)
}
