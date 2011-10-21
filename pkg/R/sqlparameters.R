##*************************************************************************
## parameters.m: set OPTIONS structure to specify default
##               parameters for sqlp.m
##
## OPTIONS.vers        : version of direction to use.  
## OPTIONS.gam         : step-length parameter,
## OPTIONS.predcorr    : whether to use Mehrotra predictor-corrector. 
## OPTIONS.expon       : exponent in decrease of centering parameter sigma. 
## OPTIONS.gaptol      : tolerance for duality gap as a fraction of the 
##                       value of the objective functions. 
## OPTIONS.inftol      : tolerance for stopping due to suspicion of 
##                       infeasibility.
## OPTIONS.steptol     : toloerance for stopping due to small steps.
## OPTIONS.maxit       : maximum number of iteration allowed 
## OPTIONS.printlevel  : 3, if want to display result in each iteration, 
##                       2, if want to display only summary,
##                       1, if want to display warning message,
##                       0, no display at all.  
## OPTIONS.scale_data  : 1, if want to scale the data before solving the problem, 
##                          else = 0
## OPTIONS.rmdepconstr : 1, if want to remove nearly dependent constraints,
##                          else = 0. 
## OPTIONS.smallblkdim : block-size threshold determining what method to compute the 
##                       schur complement matrix corresponding to semidefintie block.
##                       NOTE: this number should be small, say less than 20. 
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##*************************************************************************

sqlparameters =  function(){
  OPTIONS <- list()
  OPTIONS$vers           <- 1 
  OPTIONS$gam            <- 0  
  OPTIONS$predcorr       <- 1
  OPTIONS$expon          <- 1 
  OPTIONS$gaptol         <- 1e-8
  OPTIONS$inftol         <- 1e-8 
  OPTIONS$steptol        <- 1e-6 
  OPTIONS$maxit          <- 100
  OPTIONS$printlevel     <- 0 
  OPTIONS$scale_data     <- 0 
  OPTIONS$spdensity      <- 0.5 
  OPTIONS$rmdepconstr    <- 0 
  OPTIONS$cachesize      <- 256 
  OPTIONS$smallblkdim    <- 15
  return(OPTIONS)
}
