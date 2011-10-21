mexMatvec = function(A,y,options=0){
  m1 <- dim(A)[1]
  n1 <- dim(A)[2]
  isspA <- is(A,"sparseMatrix")
  irA <- 0
  jcA <- 0
  if(isspA){
    B <- as(A,"dgCMatrix")
    A <- B@x
    irA <- B@i
    jcA <- B@p
  }
  if(options==0)
  Z <- .C("mexMatvec_c",as.double(A),as.integer(m1),as.integer(n1),
          as.double(y),Ay=double(m1),as.integer(irA),as.integer(jcA),
          as.integer(isspA),as.integer(options))
  else
  Z <- .C("mexMatvec_c",as.double(A),as.integer(m1),as.integer(n1),
          as.double(y),Ay=double(n1),as.integer(irA),as.integer(jcA),
          as.integer(isspA),as.integer(options))
  return(Z$Ay)
}
