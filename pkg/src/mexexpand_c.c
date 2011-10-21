/***********************************************************************
* mexexpand.c : C mex file 
*
*   z = mexexpand(blk,x); 
* 
*   Input: blk   = [n1, n2, ..., nk]
*
*   Output: [x1 x1...x1  x2 x2 ... x2 .... xk xk ...xk]'    
*            n1          n2                nk
*
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
***********************************************************************/

#include <R.h>
#include <math.h>

#if !defined(MAX)
#define  MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif

/**********************************************************
* 
***********************************************************/
void mexexpand_c(int *blksize, int *numblk, double *x, double *z){    
  int k, l, blkdim, cols, idx;
  cols = 0; 
  for (k=0; k<*numblk; k++) { 
    cols = cols + blksize[k]; 
  } 
  idx = 0; 
  for (k=0; k<*numblk; k++) { 
    blkdim = blksize[k]; 
    for (l=0; l<blkdim; l++) { z[idx] = x[k]; idx++; }
  }
  return;
}
