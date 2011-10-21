%%*****************************************************************************
%% sqlp: solve an SOCP program by infeasible path-following method. 
%%
%%  [obj,X,y,Z,info,runhist] = ...
%%       sqlp(blk,At,C,b,OPTIONS,X0,y0,Z0);
%%
%%  Input: blk: a cell array describing the structure of the SOCP data.
%%          At: a cell array with At{p} 
%%         b,C: data for the SOCP instance.
%%  (X0,y0,Z0): an initial iterate (if it is not given, the default is used).
%%     OPTIONS: a structure that specifies parameters required in sqlp.m,
%%              (if it is not given, the default in sqlparameters.m is used). 
%%
%%  Output: obj  = [<C,X> <b,y>].
%%          (X,y,Z): an approximately optimal solution or a primal or dual
%%                   infeasibility certificate. 
%%          info.termcode = termination-code  
%%          info.iter     = number of iterations
%%          info.cputime  = total-time
%%          info.gap      = gap
%%          info.pinfeas  = primal_infeas
%%          info.dinfeas  = dual_infeas  
%%          runhist.pobj    = history of primal objective value. 
%%          runhist.dobj    = history of dual   objective value.
%%          runhist.gap     = history of <X,Z>. 
%%          runhist.pinfeas = history of primal infeasibility. 
%%          runhist.dinfeas = history of dual   infeasibility. 
%%          runhist.cputime = history of cputime spent.
%%          (Xiter,yiter,Ziter): last iterates.
%%
%%*************************************************************************
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%*****************************************************************************

  function [obj,X,y,Z,info,runhist] = ...
                sqlp(blk,At,C,b,OPTIONS,X0,y0,Z0);
%%                                      
%%-----------------------------------------
%% get parameters from the OPTIONS structure. 
%%-----------------------------------------
%%
   global matlabversion 
   global spdensity  iter  depconstr printlevel  solve_ok  
   global nnzmatold  matfct_options  matfct_options_old  use_LU
   global Lsymb

   warning off; 
   matlabversion = sscanf(version,'%f');
   matlabversion = matlabversion(1);
   use_LU = 0;
   
   if (nargin <= 4) | isempty(OPTIONS); OPTIONS = sqlparameters; end; 

   vers        = 1; 
   predcorr    = 1; 
   gam         = 0; 
   expon       = 1; 
   gaptol      = 1e-8;
   inftol      = 1e-8;
   steptol     = 1e-6;
   maxit       = 100;
   printlevel  = 3;
   stoplevel   = 1; 
   plotyes     = 1; 
   spdensity   = 0.5; 
   rmdepconstr = 0; 
   cachesize   = 256; 
   if isfield(OPTIONS,'vers');        vers     = OPTIONS.vers; end 
   if isfield(OPTIONS,'predcorr');    predcorr = OPTIONS.predcorr; end 
   if isfield(OPTIONS,'gam');         gam      = OPTIONS.gam; end
   if isfield(OPTIONS,'expon');       expon    = OPTIONS.expon; end
   if isfield(OPTIONS,'gaptol');      gaptol   = OPTIONS.gaptol; end
   if isfield(OPTIONS,'inftol');      inftol   = OPTIONS.inftol; end
   if isfield(OPTIONS,'steptol');     steptol  = OPTIONS.steptol; end
   if isfield(OPTIONS,'maxit');       maxit    = OPTIONS.maxit; end
   if isfield(OPTIONS,'printlevel');  printlevel  = OPTIONS.printlevel; end 
   if isfield(OPTIONS,'stoplevel');   stoplevel   = OPTIONS.stoplevel; end 
   if isfield(OPTIONS,'plotyes');     plotyes     = OPTIONS.plotyes; end 
   if isfield(OPTIONS,'spdensity');   spdensity   = OPTIONS.spdensity; end
   if isfield(OPTIONS,'rmdepconstr'); rmdepconstr = OPTIONS.rmdepconstr; end
   if isfield(OPTIONS,'cachesize');   cachesize   = OPTIONS.cachesize; end
%%
%%-----------------------------------------
%% convert matrices to cell arrays. 
%%-----------------------------------------
%%
   if ~iscell(At); At = {At}; end;
   if ~iscell(C);  C = {C}; end;
   X = X0; y = y0; Z = Z0;  
   if ~iscell(X);  X = {X}; end;
   if ~iscell(Z);  Z = {Z}; end;
%%
%%-----------------------------------------
%% validate SOCP data. 
%%-----------------------------------------
%%
   tstart = cputime; 
   [At,C,dim,numblk,X,Z] = validate(blk,At,C,b,X,y,Z);
%%
%%-----------------------------------------
%% convert unrestricted blk to linear blk. 
%%-----------------------------------------
%%
   ublkidx = zeros(size(blk,1),1); 
   for p = 1:size(blk,1) 
     if strcmp(blk{p,1},'u') & (blk{p,2} > 20)
         ublkidx(p) = 1; 
         n = 2*blk{p,2}; 
         blk{p,1} = 'l'; 
         blk{p,2} = n;
         At{p} = [At{p}; -At{p}];       
         C{p} = [C{p}; -C{p}];
         b2 = 1 + abs(b');  
         normC = 1+norm(C{p});
         normA = 1+sqrt(sum(At{p}.*At{p}));
         X{p} = max(1,max(b2./normA)) *ones(n,1);
         Z{p} = max(1,max([normA,normC])/sqrt(n)) *ones(n,1);
      end
   end
%%
%%-----------------------------------------
%% check if the matrices Ak are 
%% linearly independent. 
%%-----------------------------------------
%%
   m0 = length(b); 
   [At,b,y,indeprows,depconstr,feasible] = ...
    checkdepconstr(blk,At,b,y,rmdepconstr);
   if (~feasible)
      fprintf('\n sqlp: SOCP is not feasible'); return; 
   end
%%
%%-----------------------------------------
%% initialization
%%-----------------------------------------
%%
   A = cell(size(blk,1),1); 
   for p = 1:size(blk,1); A{p} = At{p}'; end
       
   normb = norm(b); normC = ops(C,'norm'); normA = ops(A,'norm'); 
   normX0 = ops(X0,'norm'); normZ0 = ops(Z0,'norm'); 
   m = length(b); 
   n = ops(C,'getM'); 
   AX = AXfun(blk,A,X); 
   rp = b-AX;
   Aty = Atyfun(blk,A,y);
   ZpATy = ops(Z,'+',Aty);
   ZpATynorm = ops(ZpATy,'norm');
   Rd  = ops(C,'-',ZpATy);
   obj = [blktrace(blk,C,X), b'*y];
   gap = blktrace(blk,X,Z);  
   mu  = gap/n;  
   rel_gap = gap/(1+mean(abs(obj)));
   prim_infeas = norm(rp)/(1+normb);
   dual_infeas = ops(Rd,'norm')/(1+normC);
   infeas_meas = max(prim_infeas,dual_infeas); 
   pstep = 0; dstep = 0; pred_convg_rate = 1; corr_convg_rate = 1;
   prim_infeas_bad = 0;
   termcode  = -6; 
   runhist.pobj = obj(1);
   runhist.dobj = obj(2); 
   runhist.gap  = gap;
   runhist.pinfeas = prim_infeas;
   runhist.dinfeas = dual_infeas;
   runhist.infeas  = infeas_meas;  
   runhist.step    = 0; 
   runhist.cputime = cputime-tstart; 
   ttime.preproc   = runhist.cputime; 
   ttime.pred = 0; ttime.predstep = 0; 
   ttime.corr = 0; ttime.corrstep = 0; ttime.misc = 0;
%%
%%-----------------------------------------
%% display parameters, and initial info
%%-----------------------------------------
%%
   if (printlevel>=2)
      fprintf('\n num. of constraints = %2.0d',m);      
      if dim(1); 
         fprintf('\n dim. of socp   var  = %2.0d,',dim(1)); 
         fprintf('   num. of socp blk  = %2.0d',numblk(1)); 
      end
      if dim(2); fprintf('\n dim. of linear var  = %2.0d',dim(2)); end
      if dim(3); fprintf('\n dim. of free   var  = %2.0d',dim(3)); end
      fprintf('\n********************************************');
      fprintf('***********************\n');
      fprintf('   SDPT3: Infeasible path-following algorithms'); 
      fprintf('\n********************************************');
      fprintf('***********************\n');
      if (printlevel>=3)       
         [hh,mm,ss] = mytimed(ttime.preproc);
         fprintf(' version  predcorr  gam  expon\n');
         fprintf('    NT '); 
         fprintf('     %1.0f      %4.3f',predcorr,gam);
         fprintf('   %1.0f        %1.0f\n',expon); 
         fprintf('\nit  pstep dstep p_infeas d_infeas  gap')
         fprintf('     mean(obj)    cputime\n');
         fprintf('------------------------------------------------');
         fprintf('-------------------\n');
         fprintf('%2.0f  %4.3f %4.3f %2.1e %2.1e',0,0,0,prim_infeas,dual_infeas);
         fprintf('  %2.1e %- 7.6e  %d:%d:%d',gap,mean(obj),hh,mm,ss);
      end
   end
%%
%%---------------------------------------------------------------
%% start main loop
%%---------------------------------------------------------------
%%
   [Xchol,indef(1)] = blkcholfun(blk,X); 
   [Zchol,indef(2)] = blkcholfun(blk,Z);    
   if any(indef)
      if (printlevel); fprintf('\n Stop: X or Z not positive definite'); end
      termcode = -3;
      return;
   end 
%%
   for iter = 1:maxit;  

       update_iter = 0; breakyes = 0; pred_slow = 0; corr_slow = 0; step_short = 0; 
       tstart = cputime;  
       time = zeros(1,11); 
       time(1) = cputime;
%%
%%---------------------------------------------------------------
%% predictor step.
%%---------------------------------------------------------------
%%
       if (predcorr)
          sigma = 0; 
       else 
          sigma = 1-0.9*min(pstep,dstep); 
          if (iter == 1); sigma = 0.5; end; 
       end
       sigmu = sigma*mu; 
       [par,dX,dy,dZ,coeff,L,hRd] = NTpred(blk,A,rp,Rd,sigmu,X,Z);
       if (solve_ok <= 0)
          runhist.cputime(iter+1) = cputime-tstart; 
          termcode = -4;
          break;
       end
       time(2) = cputime;
       ttime.pred = ttime.pred + time(2)-time(1);
%%
%%-----------------------------------------
%% step-lengths for predictor step
%%-----------------------------------------
%%
      if (gam == 0) 
         gamused = 0.9 + 0.09*min(pstep,dstep); 
      else
         gamused = gam;
      end 
      Xstep = steplength(blk,X,dX); 
      if (Xstep > .99e12) & (blktrace(blk,C,dX) < -1e-3) & (prim_infeas < 1e-3)
         if (printlevel); fprintf('\n Predictor: dual seems infeasible.'); end
      end
      pstep = min(1,gamused*Xstep);
      Zstep = steplength(blk,Z,dZ); 
      time(3) = cputime;        
      if (Zstep > .99e12) & (b'*dy > 1e-3) & (dual_infeas < 1e-3)
         if (printlevel); fprintf('\n Predictor: primal seems infeasible.'); end
      end
      dstep = min(1,gamused*Zstep);
      gappred = blktrace(blk,ops(X,'+',dX,pstep),ops(Z,'+',dZ,dstep)); 
      mupred  = gappred/n; 
      mupredhist(iter) = mupred; 
      ttime.predstep = ttime.predstep + time(3)-time(2);
%%
%%-----------------------------------------
%%  stopping criteria for predictor step.
%%-----------------------------------------
%%
      if (min(pstep,dstep) < steptol) & (stoplevel) 
         if (printlevel) 
            fprintf('\n  Stop: steps in predictor too short:');
            fprintf(' pstep = %3.2e,  dstep = %3.2e',pstep,dstep);
         end
         runhist.cputime(iter+1) = cputime-tstart; 
         termcode = -2; 
         breakyes = 1; 
      end
      if (iter >= 2) 
         idx = [max(2,iter-2) : iter];
         pred_slow = all(mupredhist(idx)./mupredhist(idx-1) > 0.4);
         idx = [max(2,iter-5) : iter];
         pred_convg_rate = mean(mupredhist(idx)./mupredhist(idx-1));
         pred_slow = pred_slow + (mupred/mu > 5*pred_convg_rate);
      end 
      if (~predcorr)
         if (max(mu,infeas_meas) < 1e-6) & (pred_slow) & (stoplevel)
            if (printlevel) 
               fprintf('\n  lack of progress in predictor:');
               fprintf(' mupred/mu = %3.2f, pred_convg_rate = %3.2f.',...
                         mupred/mu,pred_convg_rate);
            end
            runhist.cputime(iter+1) = cputime-tstart; 
            termcode = -1; 
            breakyes = 1;
         else 
            update_iter = 1; 
         end
      end
%%
%%---------------------------------------------------------------
%% corrector step.
%%---------------------------------------------------------------
%%
      if (predcorr) & (~breakyes)
         step_pred = min(pstep,dstep);
         if (mu > 1e-6)
            if (step_pred < 1/sqrt(3)); 
               expon_used = 1; 
            else
               expon_used = max(expon,3*step_pred^2); 
            end
         else 
            expon_used = max(1,min(expon,3*step_pred^2)); 
         end 
         sigma = min( 1, (mupred/mu)^expon_used ); 
         sigmu = sigma*mu; 
%%
         [dX,dy,dZ] = NTcorr(blk,A,par,rp,Rd,sigmu,hRd,dX,dZ,coeff,L,X,Z);
         if (solve_ok <= 0)
            runhist.cputime(iter+1) = cputime-tstart; 
            termcode = -4;
            break;
         end
         time(4) = cputime;
         ttime.corr = ttime.corr + time(4)-time(3);
%%
%%-----------------------------------
%% step-lengths for corrector step
%%-----------------------------------
%%
         if (gam == 0) 
            gamused = 0.9 + 0.09*min(pstep,dstep); 
         else
            gamused = gam;
         end            
         Xstep = steplength(blk,X,dX);
         if (Xstep > .99e12) & (blktrace(blk,C,dX) < -1e-3) & (prim_infeas < 1e-3)
            pstep = Xstep;
            if (printlevel); fprintf('\n Corrector: dual seems infeasible.'); end
         else
            pstep = min(1,gamused*Xstep);
         end
         Zstep = steplength(blk,Z,dZ);
         time(5) = cputime;
         if (Zstep > .99e12) & (b'*dy > 1e-3) & (dual_infeas < 1e-3)
            dstep = Zstep;
            if (printlevel); fprintf('\n Corrector: primal seems infeasible.'); end
         else
            dstep = min(1,gamused*Zstep);
         end   
         gapcorr = blktrace(blk,ops(X,'+',dX,pstep),ops(Z,'+',dZ,dstep)); 
         mucorr  = gapcorr/n;
         ttime.corrstep = ttime.corrstep + time(5)-time(4);
%%
%%-----------------------------------------
%%  stopping criteria for corrector step
%%-----------------------------------------
%%
         if (iter >= 2) 
            idx = [max(2,iter-2) : iter];
            corr_slow = all(runhist.gap(idx)./runhist.gap(idx-1) > 0.8); 
            idx = [max(2,iter-5) : iter];
            corr_convg_rate = mean(runhist.gap(idx)./runhist.gap(idx-1));
            corr_slow = corr_slow + (mucorr/mu > max(min(1,5*corr_convg_rate),0.8)); 
         end 
	    if (max(mu,infeas_meas) < 1e-6) & (iter > 10) & (corr_slow) & (stoplevel)
   	    if (printlevel) 
               fprintf('\n  lack of progress in corrector:');
               fprintf(' mucorr/mu = %3.2f, corr_convg_rate = %3.2f',...
                         mucorr/mu,corr_convg_rate); 
            end
            runhist.cputime(iter+1) = cputime-tstart; 
            termcode = -1; 
            breakyes = 1;
         else
            update_iter = 1;
         end
      end 
%%
%%---------------------------------------------------------------
%% udpate iterate
%%---------------------------------------------------------------
%%
      indef = [1 1]; 
      if (update_iter)
         for t = 1:10
            [Xchol,indef(1)] = blkcholfun(blk,ops(X,'+',dX,pstep)); 
            if (indef(1)); pstep = 0.8*pstep; else; break; end            
         end
	 for t = 1:10
            [Zchol,indef(2)] = blkcholfun(blk,ops(Z,'+',dZ,dstep)); 
            if (indef(2)); dstep = 0.8*dstep; else; break; end             
         end
         AXtmp = AX + pstep*AXfun(blk,A,dX);
         prim_infeasnew = norm(b-AXtmp)/(1+normb);
         if any(indef)
            if (printlevel); fprintf('\n Stop: X, Z not both positive definite'); end
            termcode = -3;
            breakyes = 1;         
         elseif (prim_infeasnew > max([rel_gap,50*prim_infeas,1e-8])) ...
            & (max(pstep,dstep)<=1) & (stoplevel)
            if (printlevel)
               fprintf('\n Stop: primal infeas deteriorated too much, %2.1e',prim_infeasnew);
            end
            termcode = -7; 
            breakyes = 1; 
         else
            X = ops(X,'+',dX,pstep);  
            y = y + dstep*dy;           
            Z = ops(Z,'+',dZ,dstep);  
         end
      end
%%---------------------------------------------------------------
%% adjust linear blk arising from unrestricted blk
%%---------------------------------------------------------------
%%
      for p = 1:size(blk,1)
         if (ublkidx(p) == 1)
            len = blk{p,2}/2;
            alpha = 0.8;
            xtmp = min(X{p}([1:len]),X{p}(len+[1:len])); 
            X{p}([1:len]) = X{p}([1:len]) - alpha*xtmp;
            X{p}(len+[1:len]) = X{p}(len+[1:len]) - alpha*xtmp;
            if (mu < 1e-4)
               Z{p} = 0.5*mu./max(1,X{p});
	    else
               ztmp = min(1,max(Z{p}([1:len]),Z{p}(len+[1:len]))); 
               beta1 = xtmp'*(Z{p}([1:len])+Z{p}(len+[1:len]));
               beta2 = (X{p}([1:len])+X{p}(len+[1:len])-2*xtmp)'*ztmp;
               beta = max(0.1,min(beta1/beta2,0.5));
               Z{p}([1:len]) = Z{p}([1:len]) + beta*ztmp;
               Z{p}(len+[1:len]) = Z{p}(len+[1:len]) + beta*ztmp;
            end
         end
      end
%%
%%---------------------------------------------------------------
%% compute rp, Rd, infeasibities, etc.
%%---------------------------------------------------------------
%%
      gap = blktrace(blk,X,Z); 
      mu  = gap/n;
      AX  = AXfun(blk,A,X); 
      rp  = b-AX;
      ZpATy = ops(Z,'+',Atyfun(blk,A,y));
      ZpATynorm = ops(ZpATy,'norm');
      Rd  = ops(C,'-',ZpATy); 
      obj = [blktrace(blk,C,X),  b'*y]; 
      rel_gap = gap/(1+mean(abs(obj))); 
      prim_infeas  = norm(rp)/(1+normb);
      dual_infeas = ops(Rd,'norm')/(1+normC);
      infeas_meas = max(prim_infeas,dual_infeas); 
      runhist.pobj(iter+1) = obj(1); 
      runhist.dobj(iter+1) = obj(2); 
      runhist.gap(iter+1)  = gap;
      runhist.pinfeas(iter+1) = prim_infeas;
      runhist.dinfeas(iter+1) = dual_infeas;
      runhist.infeas(iter+1)  = infeas_meas;
      runhist.step(iter+1)    = min(pstep,dstep); 
      runhist.cputime(iter+1) = cputime-tstart; 
      time(6) = cputime;
      ttime.misc = ttime.misc + time(6)-time(5); 
      if (printlevel>=3)
         [hh,mm,ss] = mytimed(sum(runhist.cputime));
         fprintf('\n%2.0f  %4.3f %4.3f',iter,pstep,dstep);
         fprintf(' %2.1e %2.1e  %2.1e',prim_infeas,dual_infeas,gap);
         fprintf(' %- 7.6e  %d:%d:%d',mean(obj),hh,mm,ss);
      end
%%
%%--------------------------------------------------
%% check convergence.
%%--------------------------------------------------
%%
      if (ops(X,'norm') > 1e15*normX0 | ops(Z,'norm') > 1e15*normZ0)
         termcode = 3;
         breakyes = 1; 
      end
      if (obj(2) > ZpATynorm / max(inftol,1e-13))
         termcode = 1;
         breakyes = 1;
      end
      if (-obj(1) > norm(AX) / max(inftol,1e-13))
         termcode = 2;
         breakyes = 1;
      end
      if (max(rel_gap,infeas_meas) < gaptol)
         if (printlevel)
            fprintf('\n Stop: max(relative gap, infeasibilities) < %3.2e',gaptol);
         end
         termcode = 0;
         breakyes = 1;
      end
      if (stoplevel)
         min_prim_infeas = min(runhist.pinfeas(1:iter)); 
         prim_infeas_bad = prim_infeas_bad + ...
              (prim_infeas > max(1e-10,min_prim_infeas) & (min_prim_infeas < 1e-2));
         if (mu < 1e-8)
            idx = [max(1,iter-1): iter];
         elseif (mu < 1e-4);
            idx = [max(1,iter-2): iter]; 
         else
            idx = [max(1,iter-3): iter];
         end
         idx2 = [max(1,iter-4): iter]; 
         gap_ratio2 = runhist.gap(idx2+1)./runhist.gap(idx2);
         gap_slowrate = min(0.8,max(0.6,2*mean(gap_ratio2)));
         gap_ratio = runhist.gap(idx+1)./runhist.gap(idx); 
         if (infeas_meas < 1e-4 | prim_infeas_bad) & (rel_gap < 5e-3); 
            gap_slow = all(gap_ratio > gap_slowrate) & (rel_gap < 5e-3);
            if (vers==1) 
  	       tmptol = max(prim_infeas,1e-2*dual_infeas); 
            else
  	       tmptol = max(0.5*prim_infeas,1e-2*dual_infeas); 
            end
            if (rel_gap < tmptol) 
               if (printlevel); fprintf('\n Stop: relative gap < infeasibility.'); end
               termcode = 0;
               breakyes = 1;           
            elseif (gap_slow) 
               if (printlevel); fprintf('\n Stop: progress is too slow.'); end
               termcode = -5; 
               breakyes = 1;
            end  
         elseif (prim_infeas_bad) & (iter >50) & all(gap_ratio > gap_slowrate)
            if (printlevel); fprintf('\n Stop: progress is bad.'); end
            termcode = -5;
            breakyes = 1; 
         elseif (infeas_meas < 1e-8) & (gap > 1.2*mean(runhist.gap(idx)))
            if (printlevel); fprintf('\n Stop: progress is bad.'); end
            termcode = -5;
            breakyes = 1;  
         end
         if (max(runhist.infeas) > 1e-4) & (min(runhist.infeas) < 1e-4 | prim_infeas_bad) 
            rel_gap2 = abs(diff(obj))/(1+mean(abs(obj))); 
            if (rel_gap2 < 1e-3); 
               step_short = all(runhist.step([iter:iter+1]) < 0.1) ;
            elseif (rel_gap2 < 1) 
               idx = [max(1,iter-3): iter+1];
               step_short = all(runhist.step(idx) < 0.05); 
            end
            if (step_short) 
               if (printlevel); fprintf('\n Stop: steps too short consecutively'); end
               termcode = -5; 
               breakyes = 1;      
            end
         end
      end
      if (breakyes); break; end
   end
%%---------------------------------------------------------------
%% end of main loop
%%---------------------------------------------------------------
%%
   if (termcode == -6) & (printlevel)
      fprintf('\n Stop: maximum number of iterations reached.');
   end
%%
%%---------------------------------------------------------------
%% produce infeasibility certificates if appropriate
%%---------------------------------------------------------------
%%
   if (iter >= 1) 
      param.obj         = obj;
      param.rel_gap     = rel_gap; 
      param.prim_infeas = prim_infeas;
      param.dual_infeas = dual_infeas;
      param.inftol      = inftol;
      param.m0          = m0;
      param.indeprows   = indeprows;
      param.termcode    = termcode;
      param.AX          = AX; 
      param.normX0      = normX0; 
      param.normZ0      = normZ0; 
      [X,y,Z,resid,reldist] = sqlpmisc(blk,A,C,b,X,y,Z,param); 
   end
%%
%%---------------------------------------------------------------
%% recover unrestricted blk from linear blk
%%---------------------------------------------------------------
%% 
   for p = 1:size(blk,1)
      if (ublkidx(p) == 1)
         n = blk{p,2}/2; 
         X{p} = X{p}(1:n)-X{p}(n+[1:n]); 
         Z{p} = Z{p}(1:n); 
      end
   end
%%
%%---------------------------------------------------------------
%%  print summary
%%---------------------------------------------------------------
%%
   info.termcode = termcode;
   info.iter     = iter; 
   info.gap      = gap; 
   info.pinfeas  = prim_infeas;
   info.dinfeas  = dual_infeas;
   info.cputime  = sum(runhist.cputime); 
   if (termcode == 1) | (termcode == 2)
      info.resid = resid;
   end
   nnorm.b = normb; nnorm.C = normC; nnorm.A = normA; 
   nnorm.X = ops(X,'norm'); nnorm.y = norm(y); nnorm.Z = ops(Z,'norm'); 
   sqlpsummary(runhist,ttime,termcode,resid,reldist,nnorm);
%%*****************************************************************************
