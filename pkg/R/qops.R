##********************************************************
## qops:  Fu = qops(pblk,w,f,options,u);
## 
## options = 1, Fu(i) = <wi,fi> 
##         = 2, Fu(i) = 2*wi(1)*fi(1)-<wi,fi>
##         = 3, Fui   = w(i)*fi
##         = 4, Fui   = w(i)*fi, Fui(1) = -Fui(1). 
## options = 5, Fu = w [ f'*u ; ub + fb*alp ], where 
##              alp = (f'*u + u0)/(1+f0);  
## options = 6, compute Finv*u. 
##
## Note F = w [f0  fb'; fb  I+ fb*fb'/(1+f0) ], where
##      f0*f0 - fb*fb' = 1. 
##      Finv = (1/w) [f0  -fb'; -fb  I+ fb*fb'/(1+f0) ].
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##********************************************************

qops=function(pblk,w,f,options,u){
  
  if (options >= 1 & options <= 4)
    Fu <- mexqops(pblk,w,f,options) 
  else if (options == 5){
    s <- 1 + c(0,cumsum(pblk$size)) 
    idx1 <- s[1:length(pblk$size)]
    inprod <- mexqops(pblk,f,u,1) 
    tmp <- (u[idx1]+inprod)/(1+f[idx1]) 
    Fu <- u + mexqops(pblk,tmp,f,3)
    Fu[idx1] <- inprod
    Fu <- mexqops(pblk,w,Fu,3) 
  }else if (options == 6){
    s <- 1 + c(0,cumsum(pblk$size))
    idx1 <- s[1:length(pblk$size)]
    gamprod <- mexqops(pblk,f,u,2)
    tmp <- (u[idx1]+gamprod)/(1+f[idx1]) 
    Fu <- u - mexqops(pblk,tmp,f,3) 
    Fu[idx1] <- gamprod
    Fu <- mexqops(pblk,1/w,Fu,3) 
  }
  return(Fu)
}
