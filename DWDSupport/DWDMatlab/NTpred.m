%%**********************************************************************
%% NTpred: Compute (dX,dy,dZ) for NT direction. 
%%                       
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%**********************************************************************

 function [par,dX,dy,dZ,coeff,L,hRd] = ...
          NTpred(blk,A,rp,Rd,sigmu,X,Z);

%%
%% compute NT scaling matrix
%%
    [par.gamx,par.gamz,par.dd,par.ee,par.ff] = NTscaling(blk,X,Z);
%%
%% compute schur matrix
%%
    m = length(rp); 
    schur = sparse(m,m); 
    UU = []; EE = []; Afree = []; 
    dX = cell(size(blk,1),1); dy = []; dZ = cell(size(blk,1),1); 
%%
    for p = 1:size(blk,1)
       pblk = blk(p,:); 
       if strcmp(pblk{1},'q');       
          [schur,UU,EE] = schurmat_qblk(blk,A,par,schur,UU,EE,p);
       elseif strcmp(pblk{1},'l')
          [schur,UU,EE] = schurmat_lblk(blk,A,par,schur,UU,EE,p);
       elseif strcmp(pblk{1},'u')            
          Afree = [Afree, A{p}];
       end
    end
%%
%% compute rhs
%%
    [rhs,EinvRc,hRd] = NTrhsfun(blk,A,par,X,Z,rp,Rd,sigmu);
%%
%% solve linear system
%%
    [xx,coeff,L] = linsysolve(schur,UU,Afree,EE,rhs); 
%%
%% compute (dX,dZ)
%%
    [dX,dy,dZ] = NTdirfun(blk,A,par,Rd,EinvRc,xx); 
%%**********************************************************************
