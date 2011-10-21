%%*******************************************************************
%% NTdirfun: compute (dX,dZ), given dy, for the NT direction.
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%*******************************************************************

    function [dX,dy,dZ] = NTdirfun(blk,A,par,Rd,EinvRc,xx); 
    
    global solve_ok

    dX = cell(size(blk,1),1); dZ = cell(size(blk,1),1); dy = [];
    if (any(isnan(xx)) | any(isinf(xx)))
       solve_ok = 0;
       fprintf('\n  linsysolve: solution contains NaN or inf.');
       return;
    end
%%
    m = size(A{1},1); 
    dy = xx(1:m); 
    count = m; 
%%
    for p=1:size(blk,1)
       pblk = blk(p,:);  
       if strcmp(pblk{1},'l')
          dZ{p} = Rd{p} - mexMatvec(A{p},dy,1);     
          tmp   = par.dd{p}.*dZ{p};
          dX{p} = EinvRc{p} - tmp;
       elseif strcmp(pblk{1},'q')
          dZ{p} = Rd{p} - mexMatvec(A{p},dy,1);  
          tmp = par.dd{p}.*dZ{p} + qops(pblk,qops(pblk,dZ{p},par.ee{p},1),par.ee{p},3); 
          dX{p} = EinvRc{p} - tmp;       
       elseif strcmp(pblk{1},'u'); 
          n = sum(pblk{2}); 
          dZ{p} = zeros(n,1); 
          dX{p} = xx(count+[1:n]); 
          count = count + n;  
       end
    end 
%%*******************************************************************
