%%**********************************************************************
%% blktrace: compute <X1,Z1> + ... <Xp,Zp>
%%              
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%**********************************************************************

  function trXZ = blktrace(blk,X,Z); 

  trXZ = 0; 
  for p = 1:size(blk,1)
     trXZ = trXZ + X{p}'*Z{p}; 
  end
%%**********************************************************************




