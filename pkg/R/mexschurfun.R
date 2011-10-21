mexschurfun = function(X,dd,options){
  n <- length(dd)
  isspX <- is(X,"sparseMatrix")
  irX <- 0
  jcX <- 0
  if(isspX){
    Y <- as(X,"dgTMatrix")
    irX <- Y@i
    jcX <- Y@j
    x <- Y@x
    m <- length(X)
  }else{
    m <- dim(X)[1]*dim(X)[2]
  }
  Z <- .C("mexschurfun_c",X=as.double(X),as.integer(irX),
          as.integer(jcX),as.double(dd),as.integer(n),
          as.integer(isspX),as.integer(options))
  return(Z$X)
}
