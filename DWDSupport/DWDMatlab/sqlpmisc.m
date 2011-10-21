%%*****************************************************************************
%% misc: 
%% unscale and produce infeasibility certificates if appropriate
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%*****************************************************************************

   function [X,y,Z,resid,reldist] = misc(blk,A,C,b,X,y,Z,param);

   obj         = param.obj;
   rel_gap     = param.rel_gap; 
   prim_infeas = param.prim_infeas;
   dual_infeas = param.dual_infeas;
   inftol      = param.inftol;
   m0          = param.m0;    
   indeprows   = param.indeprows;
   termcode    = param.termcode;
   AX          = param.AX;
   normX0      = param.normX0;
   normZ0      = param.normZ0;
%%
   infeas_meas = max(prim_infeas,dual_infeas);
   resid = []; reldist = [];
   Anorm = ops(A,'norm'); xnorm = ops(X,'norm'); ynorm = norm(y);
   ZpATy = ops(Z,'+',Atyfun(blk,A,y));
   ZpATynorm = ops(ZpATy,'norm');
%%
   if (termcode <= 0)
      err = min(inftol,max(infeas_meas,rel_gap));
      iflag = 0;
      if (obj(2) > 0),
         homRd = ZpATynorm/obj(2);
         if (homRd < err)
            iflag = 1;
            fprintf('\n pri_inf,dual_inf,rel_gap = %3.2e, %3.2e, %3.2e',...
                      prim_infeas,dual_infeas,rel_gap);
            termcode = 1;
         end
      elseif (obj(1) < 0),
         homrp = norm(AX)/(-obj(1)); 
         if (homrp < err) 
            fprintf('\n pri_inf,dual_inf,rel_gap = %3.2e, %3.2e, %3.2e',...
                      prim_infeas,dual_infeas,rel_gap);
            iflag = 1; 
            termcode = 2;
         end
      end
   end
   if (termcode == 1)
      fprintf('\n Stop: primal problem is suspected of being infeasible');
      rby = 1/(b'*y); y = rby*y; Z = ops(Z,'*',rby);
      resid = ZpATynorm * rby;
      reldist = ZpATynorm/(Anorm*ynorm);
   end  
   if (termcode == 2)
      fprintf('\n Stop: dual problem is suspected of being infeasible');
      tCX = blktrace(blk,C,X);
      X = ops(X,'*',1/(-tCX));
      resid = norm(AX)/(-tCX);
      reldist = norm(AX)/(Anorm*xnorm);
   end
   if (termcode == 3)
      maxblowup = max(ops(X,'norm')/normX0,ops(Z,'norm')/normZ0);
      fprintf('\n Stop: primal or dual is diverging, %3.1e',maxblowup);
   end
   if ~isempty(indeprows)
      ytmp = zeros(m0,1); 
      ytmp(indeprows) = y;
      y = ytmp; 
   end
%%*****************************************************************************
