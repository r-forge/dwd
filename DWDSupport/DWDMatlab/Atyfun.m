%%*********************************************************
%% Atyfun: compute At*y
%%
%%  Q = Atyfun(blk,A,y);
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%**********************************************************

  function Q = Atyfun(blk,A,y);

     Q = cell(size(blk,1),1);
     for p = 1:size(blk,1)
        Q{p} = (y'*A{p})'; 
     end
%%********************************************************* 

