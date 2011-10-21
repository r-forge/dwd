%%*****************************************************************************
%% checkconst: compute AAt to determine if the 
%%             constraint matrices Ak are linearly independent. 
%%              
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%*****************************************************************************

   function  [At,b,y,idxB,ndepconstr,feasible] = checkdepconstr(blk,At,b,y,rmdepconstr);
   
   global spdensity 
%%  
%% compute AAt
%%
   m = length(b);  
   AAt = sparse(m,m);  numdencol = 0; UU = []; 
   for p = 1:size(blk,1)
      pblk = blk(p,:); 
      decolidx = checkdense(At{p}'); 
      if ~isempty(decolidx); 
         n2 = size(At{p},1); 
         dd = ones(n2,1); 
         len= length(decolidx); 
         dd(decolidx) = zeros(len,1);  
         AAt = AAt + (At{p}' *spdiags(dd,0,n2,n2)) *At{p};
         tmp = At{p}(decolidx,:)'; 
         UU = [UU tmp]; 
         numdencol = numdencol + len; 
      else
         AAt = AAt + At{p}'*At{p};  
      end
   end
   if (numdencol > 0)
      fprintf('\n number of dense column in AAt = %d',numdencol); 
   end
   numdencol = size(UU,2); 
%%
%%
%%
   if ~issparse(AAt); AAt = sparse(AAt); end
   pertdiag = 1e-16*norm(AAt,'fro')*ones(m,1);
   AAt = AAt + spdiags(pertdiag,0,m,m);
   [R,indef,perm] = chol(AAt,'vector');
   if (indef)
      fprintf(' AAt is not positive definite '); 
      return;
   else
      dd = diag(R);       
   end
%%
%% find independent rows of A
%%
   idxB = [1:m]';   
   feasible = 1; ndepconstr = 0; 
   idxN = find(dd < 1e-8);
   if ~isempty(idxN)     
      fprintf('\n number of nearly dependent constraints = %1.0d',length(idxN)); 
      ndepconstr = 1; 
      if (numdencol==0)
         if (rmdepconstr)
            idxB = setdiff([1:m]',idxN);     
            fprintf('\n checkdepconstr: removing dependent constraints...');
            W = findcoeff(blk,At,idxB,idxN);
            tmp = W'*b(idxB) - b(idxN);
            nnorm = norm(tmp)/max(1,norm(b)); 
   	    tol = 1e-8;
            if (nnorm > tol)
               feasible = 0; 
               fprintf('\n checkdepconstr: inconsistent constraints exist,');
               fprintf(' problem is infeasible.');
            else
               feasible = 1; 
               for p = 1:size(blk,1) 
                  At{p} = At{p}(:,idxB);
               end
               b = b(idxB);
               y = y(idxB); 
            end
         else
            fprintf('\n To remove these constraints,');
            fprintf(' re-run sqlp.m with OPTIONS.rmdepconstr = 1.'); 
         end
      else
         fprintf('\n warning: the sparse part of AAt may be nearly');
         fprintf(' singular.');
      end
   end
%%
%%*****************************************************************************
%% findcoeff: 
%%
%% W = findcoeff(blk,At,idXB,idXN);
%% 
%% idXB = indices of independent columns of At. 
%% idxN = indices of   dependent columns of At.
%% 
%% AB = At(:,idxB); AN = At(:,idxN) = AB*W
%%
%% 
%% SDPT3: version 3.0
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%*****************************************************************************

   function W = findcoeff(blk,At,idxB,idxN);

   tol = 1e-8; 
   AB = []; AN = [];
   for p = 1:size(blk,1) 
       AB = [AB; At{p}(:,idxB)];
       AN = [AN; At{p}(:,idxN)];
   end
   [m,n] = size(AB); 
%%
%%-----------------------------------------
%% find W so that AN = AB*W
%%-----------------------------------------
%% 
   [L,U,P,Q] = lu(sparse(AB));    
   rhs  = P*AN;
   Lhat = L(1:n,:); 
   W = Q*( U \ (Lhat \ rhs(1:n,:))); 
   
   nnorm = norm(AN-AB*W,'fro')/max(1,norm(AN,'fro'));
   if (nnorm > tol) 
      fprintf('\n findcoeff: basis rows may be identified incorrectly.'); 
   end
%%*****************************************************************************
