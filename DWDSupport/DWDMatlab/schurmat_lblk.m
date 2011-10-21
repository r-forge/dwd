%%*******************************************************************
%% schurmat_lblk: 
%%
%% [schur,UU,EE] = schurmat_lblk(blk,A,par,schur,UU,EE,p);
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%*******************************************************************

   function [schur,UU,EE] = schurmat_lblk(blk,A,par,schur,UU,EE,p);
   
   n = sum(blk{p,2});  
   idxdenAl = checkdense(A{p}); 

   ddsch = par.dd{p}; 
   if ~isempty(idxdenAl); 
      idxden = idxdenAl; 
      len = length(idxden); 
      Ad = A{p}(:,idxden)*spdiags(sqrt(ddsch(idxden)),0,len,len); 
      UU = [UU, Ad];
      if isempty(EE)
         count = 0; 
      else
         count = max(max(EE(:,1)),max(EE(:,2))); 
      end
      tmp = count + [1:len]'; 
      EE = [EE; [tmp, tmp, -ones(len,1)] ]; 
      ddsch(idxden) = zeros(len,1); 
   end
   schurtmp = A{p} *spdiags(ddsch,0,n,n) *A{p}'; 
   schur = schur + schurtmp;
%%*******************************************************************
