%%*****************************************************************************
%% summary: print summary
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%*****************************************************************************

  function  summary(runhist,ttime,termcode,resid,reldist,nnorm) 

   global printlevel

   iter = length(runhist.pobj)-1; 
   obj  = [runhist.pobj(iter+1);  runhist.dobj(iter+1)];
   gap  = runhist.gap(iter+1);  
   rel_gap = gap/max(1,mean(abs(obj)));   
   prim_infeas = runhist.pinfeas(iter+1);
   dual_infeas = runhist.dinfeas(iter+1);
%%
   preproctime = ttime.preproc; 
   predtime  = ttime.pred; predsteptime = ttime.predstep; 
   corrtime  = ttime.corr; corrsteptime = ttime.corrstep; 
   misctime  = ttime.misc; 
%%
   if (printlevel >= 2) 
      fprintf('\n----------------------------------------------------\n');
      fprintf(' number of iterations   = %2.0f\n',iter);
   end
   totaltime = sum(runhist.cputime);
   if (termcode <= 0)
      if (printlevel >=2)
         fprintf(' primal objective value = %- 9.8e\n',obj(1));
         fprintf(' dual   objective value = %- 9.8e\n',obj(2));
         fprintf(' gap := trace(XZ)       = %3.2e\n',gap);
         fprintf(' relative gap           = %3.2e\n',rel_gap);
         fprintf(' actual relative gap    = %3.2e\n',-diff(obj)/(1+mean(abs(obj))));
         fprintf(' rel. primal infeas     = %3.2e\n',prim_infeas);
         fprintf(' rel. dual   infeas     = %3.2e\n',dual_infeas);
         fprintf(' norm(X), norm(y), norm(Z) = %3.1e, %3.1e, %3.1e\n',...
                   nnorm.X,nnorm.y,nnorm.Z);
         fprintf(' norm(A), norm(b), norm(C) = %3.1e, %3.1e, %3.1e\n',...
                   nnorm.A,nnorm.b,nnorm.C);
      end
   elseif (termcode == 1)
      if (printlevel >=2)
         fprintf(' residual of primal infeasibility      \n')
         fprintf(' certificate (y,Z)      = %3.2e\n',resid);
         fprintf(' reldist to infeas.    <= %3.2e\n',reldist);
      end
   elseif (termcode == 2)
      if (printlevel >=2)
         fprintf(' residual of dual infeasibility        \n')
         fprintf(' certificate X          = %3.2e\n',resid);
         fprintf(' reldist to infeas.    <= %3.2e\n',reldist);
      end
   end
   if (printlevel >=2)
      fprintf(' Total CPU time (secs)  = %3.1f  \n',totaltime);
      fprintf(' CPU time per iteration = %3.1f  \n',totaltime/iter);
      fprintf(' termination code       = %2.0f\n',termcode);
      fprintf('------------------------------------------------');
      fprintf('-------------------\n');
      fprintf(' Percentage of CPU time spent in various parts \n'); 
      fprintf('------------------------------------------------');
      fprintf('-------------------\n');
      fprintf(' preproc pred predstep corr corrstep misc\n')
      tt = [preproctime predtime predsteptime corrtime corrsteptime misctime];
      tt = tt/sum(tt)*100; 
      fprintf('   %2.1f   %2.1f    %2.1f   %2.1f   %2.1f     %2.1f\n',...
               tt(1),tt(2),tt(3),tt(4),tt(5),tt(6)); 
      fprintf('------------------------------------------------');
      fprintf('-------------------\n');
   end
%%*****************************************************************************
