%%*******************************************************************
%% NTrhsfun: compute the right-hand side vector of the 
%%           Schur complement equation for the NT direction. 
%% 
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%*******************************************************************

    function [rhs,EinvRc,hRd] = NTrhsfun(blk,A,par,X,Z,rp,Rd,sigmu,hRd,dX,dZ); 
    
    global spdensity 

    m = length(rp);   
    if (nargin > 8) 
       corrector = 1; 
    else 
       corrector = 0; 
       hRd = zeros(m,1); 
    end       
    hEinvRc = zeros(m,1); 
    rhsfree = []; 
%%
    for p = 1:size(blk,1)
       pblk = blk(p,:); 
       n = sum(pblk{2});  numblk = length(pblk{2});  
       if strcmp(pblk{1},'l')
          if (corrector)
             Rq = dX{p}.*dZ{p}; 
          else
             Rq = sparse(n,1); 
             tmp  = par.dd{p}.*Rd{p};
             tmp2 = mexMatvec(A{p},tmp); 
             hRd = hRd + tmp2;
          end
          EinvRc{p} = sigmu./Z{p}-X{p} -Rq./Z{p};
          tmp2 = mexMatvec(A{p},EinvRc{p});  
          hEinvRc = hEinvRc + tmp2;
       elseif strcmp(pblk{1},'q') 
          w = sqrt(par.gamz{p}./par.gamx{p}); 
   	  if (corrector)
             hdx = qops(pblk,w,par.ff{p},5,dX{p}); 
             hdz = qops(pblk,w,par.ff{p},6,dZ{p}); 
             hdxdz = Arrow(pblk,hdx,hdz);
             vv = qops(pblk,w,par.ff{p},5,X{p}); 
             Vihdxdz = Arrow(pblk,vv,hdxdz,1); 
             Rq = qops(pblk,w,par.ff{p},6,Vihdxdz); 
          else
             Rq = sparse(n,1); 
             tmp  = par.dd{p}.*Rd{p} + qops(pblk,qops(pblk,Rd{p},par.ee{p},1),par.ee{p},3);
             tmp2 = mexMatvec(A{p},tmp); 
             hRd = hRd + tmp2;
          end
          EinvRc{p} = qops(pblk,-sigmu./(par.gamz{p}.*par.gamz{p}),Z{p},4)-X{p}-Rq;
          tmp2 = mexMatvec(A{p},EinvRc{p});
          hEinvRc = hEinvRc + tmp2;
       elseif strcmp(pblk{1},'u') 
          rhsfree = [rhsfree; Rd{p}]; 
       end 
   end
%% 
    rhs = rp + hRd - hEinvRc; 
    rhs = [rhs; rhsfree];  
%%*******************************************************************
