%%***************************************************************
%% linsysolve: solve linear system to get dy, and direction
%%             corresponding to unrestricted variables. 
%%
%% [xx,coeff,L,M] = linsysolve(schur,UU,Afree,EE,rhs); 
%%
%% child functions: symqmr.m, linsysolvefun.m
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%***************************************************************
   
   function [xx,coeff,L] = linsysolve(schur,UU,Afree,EE,rhs); 
   
    global matlabversion 
    global spdensity  iter  depconstr  printlevel  solve_ok  
    global nnzmatold  matfct_options  matfct_options_old  use_LU
    global Lsymb

    nnzmat = 0;
    if (iter==1); use_LU = 0; end
    if isempty(nnzmatold); nnzmatold = 0; end
%%
%% diagonal perturbation, diagonally scale schur
%%
    m = length(schur); 
    diagschur = full(abs(diag(schur)));
    if (length(Afree)==0)
       if (depconstr) 
          pertdiag = 1e-15*max(1,diagschur); 
       else
          pertdiag = 1e-15*max(1e-4,diagschur); 
       end
       mexschurfun(schur,pertdiag); 
    end
%%
%% assemble coefficient matrix
%% 
    len = size(Afree,2);
    if ~isempty(EE)
       EE(:,[1 2]) = len + EE(:,[1 2]); %% adjust for ublk
    end
    EE = [(1:len)' (1:len)' zeros(len,1); EE]; 
    if isempty(EE)
       coeff.mat22 = []; 
    else
       coeff.mat22 = spconvert(EE);
    end
    coeff.mat12 = [Afree, UU]; 
    coeff.mat11 = schur; %% important to use perturbed schur matrix
    ncolU = size(coeff.mat12,2); 
%%
%% pad rhs with zero vector
%% decide which solution methods to use
%%
    rhs = [rhs; zeros(m+ncolU-length(rhs),1)]; 
    if (ncolU > 300); use_LU = 1; end
%%
%% Cholesky factorization
%%
    L = []; resnrm = []; xx = inf*ones(m,1);
    if (~use_LU)
       nnzmat = mexnnz(coeff.mat11);
       nnzmatdiff = (nnzmat ~= nnzmatold);   
       solve_ok = 1;  solvesys = 1; 
       if (nnzmat > spdensity*m^2) | (m < 500)  
          matfct_options = 'chol';
          if issparse(schur); schur = full(schur); end;
       else
           matfct_options = 'spchol';
          if ~issparse(schur); schur = sparse(schur); end;
       end
       if (printlevel); fprintf(' %s',matfct_options); end 
       if strcmp(matfct_options,'chol')
          L.matfct_options = 'chol';    
          L.perm = [1:m];
          [L.L,indef] = chol(schur); 
          if (indef)
 	         solve_ok = -2; solvesys = 0;
             fprintf('\n  chol: Schur complement matrix not pos. def.');
          end
       elseif strcmp(matfct_options,'spchol')
          if ~issparse(schur); schur = sparse(schur); end;
          L.matfct_options = 'spchol';    
          [L.L,indef,L.perm] = chol(schur,'lower','vector');
          L.Lt = L.L'; 
          if (indef)
 	         solve_ok = -2; solvesys = 0;
             fprintf('\n  chol: Schur complement matrix not pos. def.');
          end
       end    
       if (solvesys)
          if (ncolU)
             tmp = coeff.mat12'*linsysolvefun(L,coeff.mat12)-coeff.mat22; 
	         if issparse(tmp); tmp = full(tmp); end
             [L.Ml,L.Mu,L.Mp] = lu(tmp);
             tol = 1e-16; 
             idx = find(abs(diag(L.Mu)) < tol);
             if ~isempty(idx)
                pertdiag = zeros(ncolU,1); pertdiag(idx) = tol; 
                L.Mu = L.Mu + spdiags(pertdiag,0,ncolU,ncolU); 
		        fprintf('*');
             end
          end
          [xx,resnrm,solve_ok] = symqmr(coeff,rhs,L);
          if (solve_ok<=0) & (printlevel)
             fprintf('\n  warning: symqmr fails: %3.1f.',solve_ok); 
          end
          if (printlevel>=3); fprintf(' %2.0d',length(resnrm)-1); end
       end
       if (solve_ok < 0) | (solvesys == 0)
          if (m < 5000) 
             use_LU = 1;
             if (printlevel)
                fprintf('\n  switch to LU factor.'); 
             end
          elseif (solve_ok==-0.5)
             solve_ok = 1; 
          end
       end
    end
%%
%% symmetric indefinite or LU factorization
%%
    if (use_LU)
       nnzmat = mexnnz(coeff.mat11)+mexnnz(coeff.mat12); 
       nnzmatdiff = (nnzmat ~= nnzmatold);  
       solve_ok = 1; 
       if ~isempty(coeff.mat22)
          raugmat = [coeff.mat11, coeff.mat12; coeff.mat12', coeff.mat22]; 
       else
          raugmat = coeff.mat11; 
       end
       if (nnzmat > spdensity*m^2) | (m+ncolU < 500) 
          matfct_options = 'lu';     
       else
          matfct_options = 'splu';
       end
       if (printlevel); fprintf(' %s ',matfct_options); end 
       if strcmp(matfct_options,'lu') 
          if issparse(raugmat); raugmat = full(raugmat); end
          [L.l,L.u,L.p] = lu(raugmat); 
          L.matfct_options = 'lu'; 
          L.p = sparse(L.p); 
          idx = find(abs(diag(L.u)) < 1e-15); 
          if ~isempty(idx)
             if (printlevel); fprintf('\n  lu: modify diag(L.u)'); end
             n = length(raugmat); lam = zeros(n,1); lam(idx) = 1e-15; 
             L.u = L.u + spdiags(lam,0,n,n);
          end
          [ii,jj] = find(L.p); [dummy,idx] = sort(ii); L.perm = jj(idx); 
       end
       if strcmp(matfct_options,'splu') 
          if ~issparse(raugmat); raugmat = sparse(raugmat); end  
          L.perm = [1:length(raugmat)];  
          L.matfct_options = 'splu';  
          %try,
          %  [L.l,L.u,L.p,L.q] = lu(raugmat);
          %  [ii,jj] = find(L.q); L.q = ii;
          %catch
          %  L.q = colamd(raugmat);
          %  [L.l,L.u,L.p] = lu(raugmat(:,L.q));
          %end
          %[ii,jj] = find(L.p); L.p = ii;
          [L.l,L.u,L.p,L.q] = lu(raugmat);
       end
       [xx,resnrm,solve_ok] = mybicgstab(coeff,rhs,L);       
       if (solve_ok<=0) & (printlevel)
          fprintf('\n  warning: mybicgstab fails: %3.1f,',solve_ok); 
       end
       if (printlevel>=3); fprintf(' %2.0d',length(resnrm)-1); end
    end
%%
    nnzmatold = nnzmat; 
    matfct_options_old = matfct_options; 
%%***************************************************************
