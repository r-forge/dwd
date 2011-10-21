%%***************************************************************************
%% steplength: compute xstep such that  X + xstep*dX >= 0.
%%
%% [xstep] = steplength(blk,X,dX);
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%***************************************************************************

    function xstep = steplength(blk,X,dX);

%%
    for p = 1:size(blk,1)
       pblk = blk(p,:); 
       numblk = length(pblk{2}); 
       if (any(isnan(dX{p})) | any(isinf(dX{p}))); xstep = 0; break; end; 
       if strcmp(pblk{1},'q')
          aa = qops(pblk,dX{p},dX{p},2); 
          bb = qops(pblk,dX{p},X{p},2); 
          cc = qops(pblk,X{p},X{p},2);
          dd = bb.*bb - aa.*cc; 
          tmp = min(aa,bb); 
          idx = find(dd > 0 & tmp < 0); 
          steptmp = 1e12*ones(numblk,1); 
          if ~isempty(idx)
             steptmp(idx) = -(bb(idx)+sqrt(dd(idx)))./aa(idx);       
          end
          idx = find(abs(aa) < eps & bb < 0); 
          if ~isempty(idx)
             steptmp(idx) = -cc(idx)./(2*bb(idx)); 
          end
          %%
          %% also need first component to be non-negative
          %%
          ss = 1 + [0, cumsum(pblk{2})];
          ss = ss(1:length(pblk{2})); 
          dX0 = dX{p}(ss); 
          X0 = X{p}(ss); 
          idx = find(dX0 < 0 & X0 > 0); 
          if ~isempty(idx)
             steptmp(idx) = min(steptmp(idx),-X0(idx)./dX0(idx)); 
          end
          xstep(p) = min(steptmp); 
       elseif strcmp(pblk{1},'l')
          idx = find(dX{p} < 0); 
          if ~isempty(idx)
             xstep(p) = min(-X{p}(idx)./dX{p}(idx));  
          else 
             xstep(p) = 1e12;
          end
       elseif strcmp(pblk{1},'u')
          xstep(p) = 1e12; 
       end
    end
    xstep = min(xstep); 
%%***************************************************************************

