/**********************************************************************
*  mextriangsp: given upper triangular U,
*  options = 1 (default), solves     U *y = b (backward substitutions)
*          = 2          , solves     U'*y = b (forward substitutions). 
*
*  Important: U is assumed to be sparse. 
*
*  y = mextriangsp(Uinput,b,options)
*
*  If options = 1, Uinput = transpose of U. 
*   
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
*********************************************************************/

#include "R.h"
#include <math.h>

/*************************************************************
  ubsolve2: solve Ux = b for x by backward substitutions. 
  Ut: sparse, Ut = transpose(U)
**************************************************************/
void ubsolve2(int n, const double *Ut, const int *irUt, const int *jcUt, 
              double *b, double *x)

{ int j, k, kstart, kend, idx; 
  double tmp; 
  
  /************************************************
      x[j]*Ut[j,j] + x[j+1:n]*Ut[j+1:n,j] = b[j]
  ************************************************/

  x[n-1] = b[n-1]/Ut[jcUt[n]-1];
 
  for (j=n-2; j>=0; j--) {
      kstart = jcUt[j]+1; kend = jcUt[j+1]; 
      tmp = 0.0; 
      for (k=kstart; k<kend; k++) { 
	  idx = irUt[k]; 
          tmp += Ut[k]*x[idx];  	 
      }
      x[j] = (b[j]-tmp)/Ut[kstart-1]; 
  }
  return; 
}
/*************************************************************
  lbsolve2: solve U'*x = b for x by forward substitutions. 

**************************************************************/
void lbsolve2(int n, const double *U, const int *irU, const int *jcU, 
              double *b, double *x)

{ int j, k, kstart, kend, idx; 
  double tmp; 
  
  /************************************************
      x[j]*U[j,j] + x[0:j-1]*U[0:j-1,j] = b[j]
  ************************************************/
  x[0] = b[0]/U[0]; 
  for (j=1; j<n; j++) {
      kstart = jcU[j]; kend = jcU[j+1]-1; 
      tmp = 0.0; 
      for (k=kstart; k<kend; k++) { 
	  idx = irU[k]; 
          tmp += U[k]*x[idx];  	 
      }
      x[j] = (b[j]-tmp)/U[kend]; 
  }
  return; 
}

/*************************************************************
*   PROCEDURE mextriangsp_c - Entry for R
**************************************************************/
void mextriangsp_c(double *U, int *irU, int *jcU, double *b, double *x,
		   int *n, int *options){
  if (*options==1){ 
    ubsolve2(*n,U,irU,jcU,b,x);
  } else if (*options==2) { 
    lbsolve2(*n,U,irU,jcU,b,x);
  }   
  return;
}

