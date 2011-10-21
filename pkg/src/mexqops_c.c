/***********************************************************************
* mexqops.c : C mex file 
*
*   z = mexqops(blk,x,y,options);
*
*   Input: blk   = [n1, n2, ... nk]
*          x = n-vector or k-vector
*          y = n-vector, where n = n1+...+nk
*  
*   options = 1, z = k-vector, z(i) = <xi,yi> 
*           = 2, z = k-vector, z(i) = 2*xi(1)yi(1) - <xi,yi>
*           = 3, z = n-vector, zi   = x(i)*yi 
*           = 4, z = n-vector, zi   = x(i)*yi, zi(1) = -zi(1). 
*
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
***********************************************************************/

#include <R.h>
#include <math.h>

/**********************************************************
* ops1 
**********************************************************/
void ops1(double *x, double *y, double *z, int numblk, 
	  int *cumblk, int options) {  
  int  j, l, jstart, jend;
  double tmp; 
   
  if (options == 1) { 
    for (l=0; l<numblk; ++l) {  
      jstart = cumblk[l]; jend = cumblk[l+1];
      tmp = 0; 
      for (j=jstart; j<jend; ++j) { tmp += x[j]*y[j]; }
      z[l] = tmp; 
    }
  }
  else if (options == 2) { 
    for (l=0; l<numblk; ++l) {  
      jstart = cumblk[l]; jend = cumblk[l+1];
      tmp = x[jstart]*y[jstart]; 
      for (j=jstart+1; j<jend; ++j) { tmp -= x[j]*y[j]; }
      z[l] = tmp;
    }          
  }
  return;
}

/**********************************************************
* ops3 
**********************************************************/
void ops3(double *x, double *y, double *z, int numblk, 
	  int *cumblk, int options) {  
  int  j, l, jstart, jend;
  double tmp;
   
  if (options == 3) {
    for (l=0; l<numblk; ++l) {  
      jstart = cumblk[l]; jend = cumblk[l+1];
      tmp = x[l];
      for (j=jstart; j<jend; ++j) { 
	z[j] = tmp*y[j]; }
    }
  }
  else if (options == 4) { 
    for (l=0; l<numblk; ++l) {  
      jstart = cumblk[l]; jend = cumblk[l+1];
      tmp = x[l];
      for (j=jstart; j<jend; ++j) { 
	z[j] = tmp*y[j]; }
      z[jstart] = -z[jstart];
    }
  }
  return;
}

void mexqops_c(double *x, double *y, double *z, int *numblk, 
	     int *cumblk, int *options){
  int blk = *numblk, opt = *options;
  if(opt==1||opt==2)
    ops1(x, y, z, blk, cumblk, opt);
  else if(opt==3||opt==4)
    ops3(x, y, z, blk, cumblk, opt);
  return;
}

