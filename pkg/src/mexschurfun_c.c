/************************************************************************
*  mexschurfun(X,dd,options)         
*  options = 1, add dd to the diagonal of X (a square matrix) 
*  options = 2, diagonally scale X to diag(dd)*X*diag(dd). 
*
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
************************************************************************/

#include "R.h"
#include <math.h>

/********************************************************************
  PROCEDURE mexFunction - Entry for Matlab
*********************************************************************/
void mexschurfun_c(double *X, int *irX, int *jcX, double *dd, 
		   int *n, int *isspX, int *options){
  int    j, jn, k, kstart, kend, r; 
  double tmp, tmp2;
  
  if (*options==1) { 
    if (isspX) {
      for (j=0; j<*n; j++) {
	kstart = jcX[j]; kend = jcX[j+1]; 
	for (k=kstart; k<kend; k++) { 
	  r = irX[k];
	  if (r==j) { X[k] += dd[j]; break; } 
	}
      }
    } else { 
        for (j=0; j<*n; j++) { jn = j*(*n); X[j+jn] += dd[j]; }
      }
  } else {
     if (isspX) {
        for (j=0; j<*n; j++) {
           kstart = jcX[j]; kend = jcX[j+1]; 
           tmp = dd[j]; 
           for (k=kstart; k<kend; k++) { 
	      r = irX[k];
	      tmp2 = tmp*dd[r]*X[k];
              X[k] = tmp2; 
	   }
	}
      } else { 
        for (j=0; j<*n; j++) { 
           jn = j*(*n);
           tmp = dd[j]; 
           for (k=0; k<*n; k++) {  
	       tmp2 = tmp*dd[k]*X[k+jn];                
               X[k+jn] = tmp2;
	   }
	}
      }
  }
  return;
}

