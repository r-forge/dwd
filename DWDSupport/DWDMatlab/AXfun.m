%%*********************************************************
%% AXfun: compute A*X
%%
%%   AX = AXfun(blk,A,X);
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%**********************************************************

  function AX = AXfun(blk,A,X);

  m = size(A{1},1); 
  AX = zeros(m,1); 

  for p = 1:size(blk,1); 
     AX = AX + A{p}*X{p}; 
  end
%%*********************************************************
