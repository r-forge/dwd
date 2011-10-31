##*****************************************************************************
## misc: 
## unscale and produce infeasibility certificates if appropriate
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##*****************************************************************************

sqlpmisc=function(blk,A,At,C,b,X,y,Z,param){

  obj         <- param$obj
  rel_gap     <- param$rel_gap 
  prim_infeas <- param$prim_infeas
  dual_infeas <- param$dual_infeas
  inftol      <- param$inftol
  m0          <- param$m0    
  indeprows   <- param$indeprows
  termcode    <- param$termcode
  AX          <- param$AX
  normX0      <- param$normX0
  normZ0      <- param$normZ0

  infeas_meas <- max(prim_infeas,dual_infeas)
  resid <- NULL
  reldist <- NULL
#  nargin <<- 2
  Anorm <- ops(A,"norm")
  xnorm <- ops(X,"norm")
  ynorm <- normsvd(y)
#  nargin <<- 3
  ZpATy <- ops(Z,"+",Atyfun(blk,At,y))
#  nargin <<- 2
  ZpATynorm <- ops(ZpATy,"norm")

  if (termcode <= 0){
    err <- min(inftol,max(infeas_meas,rel_gap))
    iflag <- 0
    if (obj[2] > 0){
      homRd <- ZpATynorm/obj[2]
      if (homRd < err){
        iflag <- 1
        print(c("pri_inf,dual_inf,rel_gap =",
                prim_infeas,dual_infeas,rel_gap))
        termcode <- 1
      }
    }else if (obj[1] < 0){
      homrp <- normsvd(AX)/(-obj[1]) 
      if (homrp < err) {
        print(c(" pri_inf,dual_inf,rel_gap = ",
                prim_infeas,dual_infeas,rel_gap))
        iflag <- 1
        termcode <- 2
      }
    }
  }
  if (termcode == 1){
    print("Stop: primal problem is suspected of being infeasible")
    rby <- 1/(t(b)%*%y)
    y <- rby*y
#    nargin <<- 2
    Z <- ops(Z,"*",rby)
    resid <- ZpATynorm * rby
    reldist <- ZpATynorm/(Anorm*ynorm)
  }
  if (termcode == 2){
    print("Stop: dual problem is suspected of being infeasible")
    tCX <- blktrace(blk,C,X)
#    nargin <<- 3
    X <- ops(X,"*",1/(-tCX))
    resid <- normsvd(AX)/(-tCX)
    reldist <- normsvd(AX)/(Anorm*xnorm)
  }
  if (termcode == 3){
#    nargin <<- 2
    maxblowup <- max(ops(X,"norm")/normX0,ops(Z,"norm")/normZ0)
    print(c(" Stop: primal or dual is diverging",maxblowup))
  }
  if(!is.null(indeprows)){
    ytmp <- zeros(m0,1) 
    ytmp[indeprows] <- y
    y <- ytmp
  }
  return(list(X=X,y=y,Z=Z,resid=resid,reldist=reldist))
}
