%%*******************************************************************
%% schurmat_qblk: compute schur matrix corresponding to SOCP blocks.
%%
%% NT  direction: output = schur + Ae*Ae' - Ad*Ad'
%%
%% where schur = A*D*A', and Ad is the modification to ADA' 
%% so that the latter is positive definite. 
%%
%% [schur,UU,EE] = schurmat_qblk(blk,A,schur,UU,EE,p);
%% 
%% UU: stores the dense columns of Ax, Ae, Ad, and possibly 
%%     those of A*D^{1/2}. It has the form UU = [Ax Ae Ad]. 
%% EE: stores the assocaited (2,2) block matrix when the
%%     output matrix is expressed as an augmented matrix.
%%     It has the form EE = [0 -lam 0; -lam 0 0; 0 0 I].
%% 
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%*******************************************************************

   function [schur,UU,EE] = schurmat_qblk(blk,A,par,schur,UU,EE,p);
      
   if isempty(EE) 
      count = 0; 
   else 
      count = max(max(EE(:,2)),max(EE(:,1))); 
   end
   pblk = blk(p,:); n = sum(pblk{2}); numblk = length(pblk{2}); 
%%   
   ddsch = par.dd{p};       
   Ae = qprod(pblk,A{p},par.ee{p}); 
   idxden = checkdense(Ae);
   if ~isempty(idxden)
      spcolidx = setdiff([1:numblk],idxden); 
      s = 1 + [0 cumsum(pblk{2})];
      idx = s(idxden); 
      tmp = zeros(n,1); 
      tmp(idx) = sqrt(2*abs(ddsch(idx))); 
      Ad = qprod(pblk,A{p},tmp);              
      ddsch(idx) = abs(ddsch(idx)); 
      len = length(idxden);
      w2 = par.gamz{p}./par.gamx{p}; 
      lam = w2(idxden); 
      UU = [UU, Ae(:,idxden)*spdiags(sqrt(lam),0,len,len), Ad(:,idxden)]; 
      tmp = count+[1:len]'; 
      EE = [EE; [tmp, tmp, -lam; len+tmp, len+tmp, ones(len,1)] ]; 
      count = count + 2*len; 
      Ae = Ae(:,spcolidx);      
      schur = schur + (Ae*Ae');
   else
      schur = schur + (Ae*Ae');
   end
 
   idxdenAq = checkdense(A{p}); 
   if ~isempty(idxdenAq);
      idxden = idxdenAq;  
      len = length(idxden);               
      Ad = A{p}(:,idxden)*spdiags(sqrt(abs(ddsch(idxden))),0,len,len); 
      UU = [UU, Ad];
      tmp = count+[1:len]'; 
      EE = [EE; [tmp, tmp, -sign(ddsch(idxden))]]; 
      count = count + len; 
      ddsch(idxden) = zeros(len,1); 
   end  
   schurtmp = A{p} *spdiags(ddsch,0,n,n) *A{p}'; 
   schur = schur + schurtmp;
%%*******************************************************************
