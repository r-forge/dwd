mextriang = function(U,b,options=0){
  n <- length(b)
  Z <- .C("mextriang_c",as.double(U),as.double(b),
          y=double(n),as.integer(n),as.integer(options))
  return(Z$y)
}
