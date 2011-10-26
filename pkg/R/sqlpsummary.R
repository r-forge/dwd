##*****************************************************************************
## summary: print summary
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##*****************************************************************************

sqlpsummary = function(runhist,ttime,termcode,resid,reldist,nnorm,global_var){

  printlevel = global_var$printlevel
  iter <- length(runhist$pobj)-1; 
  obj  <- rbind(runhist$pobj[iter+1],runhist$dobj[iter+1])
  gap  <- runhist$gap[iter+1]  
  rel_gap <- gap/max(1,mean(abs(obj)))   
  prim_infeas <- runhist$pinfeas[iter+1]
  dual_infeas <- runhist$dinfeas[iter+1]

  preproctime <- ttime$preproc 
  predtime  <- ttime$pred
  predsteptime <- ttime$predstep
  corrtime  <- ttime$corr
  corrsteptime <- ttime$corrstep
  misctime  <- ttime$misc

  if (printlevel >= 2) 
    print(c(" number of iterations   = ",iter))

  totaltime <- sum(runhist$cputime)
  if (termcode <= 0){
    if (printlevel >=2){
      print(c(" primal objective value = ",obj[1]))
      print(c(" dual   objective value = ",obj[2]))
      print(c(" gap := trace(XZ)       = ",gap))
      print(c(" relative gap           = ",rel_gap))
      print(c(" actual relative gap    = ",-diff(obj)/(1+mean(abs(obj)))))
      print(c(" rel. primal infeas     = ",prim_infeas))
      print(c(" rel. dual   infeas     = ",dual_infeas))
      print(c(" norm(X), norm(y), norm(Z) = ",nnorm$X,nnorm$y,nnorm$Z))
      print(c(" norm(A), norm(b), norm(C) = ",nnorm$A,nnorm$b,nnorm$C))
    }
  }else if (termcode == 1){
    if (printlevel >=2){
      print(" residual of primal infeasibility  ")
      print(c(" certificate (y,Z)      = ",resid))
      print(c(" reldist to infeas.    <= ",reldist))
    }
  }else if (termcode == 2){
    if (printlevel >=2){
      print(" residual of dual infeasibility     ")
      print(c(" certificate X          = ",resid))
      print(c(" reldist to infeas.    <= ",reldist))
    }
  }
  if (printlevel >=2){
    print(c(" Total CPU time (secs)  = ",totaltime))
    print(c(" CPU time per iteration = ",totaltime/iter))
    print(c(" termination code       = ",termcode))
    print(" Percentage of CPU time spent in various parts")
    print(" preproc pred predstep corr corrstep misc")
    tt <- c(preproctime,predtime,predsteptime,corrtime,corrsteptime,misctime)
    tt <- tt/sum(tt)*100 
    print(c(tt[1],tt[2],tt[3],tt[4],tt[5],tt[6])) 
  }
}
