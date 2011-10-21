%%******************************************************************
%% blkcholfun: compute Cholesky factorization of X. 
%%          
%%  [Xchol,indef] = blkcholfun(blk,X); 
%%  
%%  X = Xchol'*Xchol;
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%******************************************************************
 
  function [Xchol,indef] = blkcholfun(blk,X); 

  if ~iscell(X); 
     indef = 0; 
     n = length(X); 
     if strcmp(blk{1},'q') 
        gamx = mexqops(blk{2},X,X,2); 
        if any(gamx <= 0) 
           indef = 1; 
        end
        Xchol = [];   
     elseif strcmp(blk{1},'l'); 
        if any(X <= 0) 
           indef = 1;
        end 
        Xchol = [];  
     elseif strcmp(blk{1},'u')
        Xchol = [];
     end
  else 
     Xchol = cell(size(X));  
     for p = 1:size(blk,1) 
        [Xchol{p},indef(p)] = blkcholfun(blk(p,:),X{p});
     end
     indef = max(indef); 
  end 
%%******************************************************************
