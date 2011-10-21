%%*************************************************************************
%% linsysolvefun: Solve H*x = b
%%
%% x = linsysolvefun(L,b)
%% where L contains the triangular factors of H. 
%% 
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%*************************************************************************
 
  function x = linsysolvefun(L,b)
 
  x = zeros(size(b)); 
  for k=1:size(b,2)
     if strcmp(L.matfct_options,'chol')
        x(L.perm,k) = mextriang(L.L, mextriang(L.L,b(L.perm,k),2), 1); 
     elseif strcmp(L.matfct_options,'spchol')
        x(L.perm,k) = mextriangsp(L.Lt,mextriangsp(L.L,b(L.perm,k),2), 1); 
     elseif strcmp(L.matfct_options,'lu')
        x(:,k) = L.u \ (L.l \ b(L.perm,k));
     elseif strcmp(L.matfct_options,'splu')     
        x(L.perm,k) = L.q*( L.u \ (L.l \ (L.p*b(L.perm,k))));
     end
  end
%%*************************************************************************
