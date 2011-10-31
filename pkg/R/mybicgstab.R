##*************************************************************************
## mybicgstab: bicgstab
##  
## [xx,resnrm,solve_ok] = mybicgstab(A,b,M1,tol,maxit)
##
## iterate on  bb - (M1)*AA*x
##
## r = b-A*xtrue; 
##
## SDPT3: version 3.1 
## Copyright (c) 1997 by
## K.C. Toh, M.J. Todd, R.H. Tutuncu
## Last Modified: 15 Sep 2004
##*************************************************************************

mybicgstab=function(A,b,M1,tol,maxit){

  resnrm <- NULL
  N <- length(b) 
  if (nargs() < 5){
    if(is.null(A$mat22)){
      maxit <- 5
    }else{  
      maxit <- min(50,max(5,max(dim(A$mat22))))
    }
  }
  if (nargs() < 4)
    tol = 1e-8
  tolb <- min(1e-4,tol*normsvd(b))
  flag <- 1
 
  x <- zeros(N,1) 
  r <- b - matvecsp(A,x)
  err <- normsvd(r)
  resnrm[1] <- err
  minresnrm <- err
  xx <- x  
  if (err < tolb)
    return

  omega <- 1.0
  p <- zeros(N,1)
  v <- zeros(N,1) 
  alp <- 0
  rho_1 <- 0
  r_tld <- r

  smtol <- 1e-40 
  for(it in 1:maxit){           
    rho <- t(r_tld)%*%r
    rho <- rho[1,1]
    if (abs(rho) < smtol){
      flag <- 2
      break 
    }
    if (it > 1){
      beta <- (rho/rho_1)* (alp/omega)
      p <- r + beta*(p - omega*v)
    }else{
      p <- r
    }
    p_hat <- linsysolvefun(M1,p) 
    v <- matvecsp(A,p_hat)
    alp <- rho / (t(r_tld)%*%v)
    alp <- alp[1,1]
    s <- r - alp*v
    s_hat <- linsysolvefun(M1,s)
    t <- matvecsp(A,s_hat)
    omega <- t(t)%*%s / (t(t)%*%t)
    omega <- omega[1,1]
    x <- x + alp*p_hat + omega*s_hat              
    r <- s - omega*t
    rho_1 <- rho

    ## check convergence

    err <- normsvd(r)
    resnrm[it+1] <- err
    minresnrm <- min(minresnrm,err) 
    if (err < tolb)
      break;  
    
    if ((err > 1e5*minresnrm) & (it > 50)){
      flag <- -0.5
      break
    } 
    if (abs(omega) < smtol){
      flag <- 2
      break 
    }
  }
  xx <- x
  
  return(list(xx=xx,resnrm=resnrm,solve_ok=flag))
}

##*************************************************************************
## matvecsp: matrix-vector multiply.
## matrix = [A.mat11 A.mat12; A.mat12' A.mat22]
##*************************************************************************

matvecsp = function(A,x){
  
  m <- dim(A$mat11)[1]
  m2 <- length(x)-m 
  if (m2 > 0)
    x1 = x[1:m] 
  else
    x1 = x
  Ax <- mexMatvec(A$mat11,x1,0)
  if (m2 > 0){
    x2 <- x[m+c(1:m2)]
    Ax <- Ax + mexMatvec(A$mat12,x2,0)
    Ax2 <- mexMatvec(A$mat12,x1,1) + mexMatvec(A$mat22,x2,0)
    Ax = c(Ax,Ax2)  
  }
  return(Ax)
}

