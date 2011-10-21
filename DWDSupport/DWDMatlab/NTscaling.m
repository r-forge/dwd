%%**********************************************************************
%% NTscaling: Compute NT scaling matrix
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%**********************************************************************

 function [gamx,gamz,dd,ee,ff] = NTscaling(blk,X,Z);

    numblk = size(blk,1);
    gamx = cell(numblk,1); gamz = cell(numblk,1); 
    dd = cell(numblk,1); ee = cell(numblk,1); ff = cell(numblk,1);   
%%
    for p = 1:size(blk,1)
       pblk = blk(p,:); 
       numblk = length(pblk{2});  
       n = sum(pblk{2});  
       if strcmp(pblk{1},'l')
          dd{p} = X{p}./Z{p};       
       elseif strcmp(pblk{1},'q');  
          gamx{p} = sqrt(qops(pblk,X{p},X{p},2)); 
          gamz{p} = sqrt(qops(pblk,Z{p},Z{p},2)); 
          w2 = gamz{p}./gamx{p};  w = sqrt(w2); 
          dd{p} = qops(pblk,1./w2,ones(n,1),4);
          tt = qops(pblk,1./w,Z{p},3) - qops(pblk,w,X{p},4);
          gamtt = sqrt(max(1e-16,qops(pblk,tt,tt,2))); 
          ff{p} = qops(pblk,1./gamtt,tt,3); 
          ee{p} = qops(pblk,sqrt(2)./w,ff{p},4); 
       end
    end
%%**********************************************************************




