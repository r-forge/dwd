##************************************************************************
## NTcorr: corrector step for the NT direction. 
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##************************************************************************

NTcorr = function(blk,A,par,rp,Rd,sigmu,hRd,dX,dZ,coeff,L,X,Z,global_var){

  matfct_options = global_var$matfct_options
  printlevel = global_var$printlevel
  myrhs <- NTrhsfun(blk,A,par,X,Z,rp,Rd,sigmu,hRd,dX,dZ)
  rhs <- myrhs$rhs
  EinvRc <- myrhs$EinvRc
  m <- length(rp)
  ncolU <- dim(coeff$mat12)[2] 
  if(is.null(ncolU))
    ncolU <- 0
  rhs <- rbind(rhs,zeros(m+ncolU-length(rhs),1)) 

  solve_ok <- 1
  resnrm <- normsvd(rhs)
  if((matfct_options=="chol")|(matfct_options=="spchol")){
    myqmr <-  symqmr(coeff,rhs,L)
    xx <- myqmr$xx
    resnrm <- myqmr$resnrm
    solve_ok <- myqmr$solve_ok
    #if ((solve_ok<=0) & (printlevel))
    #  print('warning: symqmr fails: ',solve_ok); 
  }else{
    mybicg <- mybicgstab(coeff,rhs,L)
    xx <- mybicg$xx 
    resnrm <- mybicg$resnrm
    solve_ok <- mybicg$solve_ok
    if ((solve_ok<=0) & (printlevel))
      print(c("warning: mybicgstab fails:",solve_ok)) 
  }
  if(printlevel>=3)
    print(length(resnrm)-1)
  if ((any(is.nan(xx)) | any(is.infinite(xx)))){
    solve_ok <- 0
    print("NTcorr: dy contains NaN or inf.")
  }
  global_var$solve_ok = solve_ok
  mydir <- NTdirfun(blk,A,par,Rd,EinvRc,xx,global_var)
  dX <- mydir$dX
  dy <- mydir$dy
  dZ <- mydir$dZ
  global_var = mydir$global_var

  return(list(dX=dX,dy=dy,dZ=dZ,global_var=global_var))
}
