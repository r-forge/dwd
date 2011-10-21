%%*************************************************************************
%% symqmr: symmetric QMR with left (symmetric) preconditioner. 
%%         The preconditioner used is based on the analytical
%%         expression of inv(A).  
%%
%% [x,resnrm,solve_ok] = symqmr(A,b,L,tol,maxit) 
%%
%% child function: linsysolvefun.m 
%%
%% A = [mat11 mat12; mat12' mat22].
%% b = rhs vector.
%% if matfct_options = 'chol' or 'spchol' 
%%    L = Cholesky factorization of (1,1) block. 
%%    M = Cholesky factorization of 
%%        Schur complement of A ( = mat12'*inv(mat11)*mat12-mat22).
%% else
%%    L = triangular factors of A.
%%    M = not relevant.
%% end
%% resnrm = norm of qmr-generated residual vector b-Ax. 
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%*************************************************************************

   function  [xx,resnrm,solve_ok] = symqmr(A,b,L,tol,maxit) 

   N = length(b); 
   if (nargin < 5); maxit = max(50,max(5,length(A.mat22))); end;
   if (nargin < 4); tol = 1e-10; end; 
   tolb = min(1e-4,tol*norm(b));

   solve_ok = 1; 
   x = zeros(N,1);
   if isstruct(A); Aq = matvec(A,x); else; Aq=mexMatvec(A,x); end;     
   r = b-Aq;  
   err = norm(r); resnrm(1) = err; minres = err; xx = x; 
   if (err < tolb); return; end         

   q = precond(A,L,r); 
   tau_old   = norm(q);      
   rho_old   = r'*q; 
   theta_old = 0; 
   d = zeros(N,1); 
   res = r; Ad = zeros(N,1);
%%      
%% main loop
%%
   tiny = 1e-30; 
   for iter = 1:maxit 

       if isstruct(A); Aq = matvec(A,q); else; Aq=mexMatvec(A,q); end;     
       sigma = q'*Aq; 
       if (abs(sigma) < tiny)
          solve_ok = 2; break;
       else
          alpha = rho_old/sigma; 
          r = r - alpha*Aq;
       end
       u = precond(A,L,r); 

       theta = norm(u)/tau_old; c = 1/sqrt(1+theta^2); 
       tau = tau_old*theta*c;
       gam = (c^2*theta_old^2); eta = (c^2*alpha); 
       d = gam*d + eta*q;
       x = x + d; 
%%
       Ad = gam*Ad + eta*Aq;
       res = res - Ad; 
       err = norm(res); resnrm(iter+1) = err; 
       if (err < minres); xx = x; minres = err; end
       if (err < tolb); break; end        
       if (iter > 5) 
          if (err > 0.98*mean(resnrm(iter-5:iter)))
             solve_ok = -0.5; break; 
          end
       end
%% 
       if (abs(rho_old) < tiny)
          solve_ok = 2; break;
       else
          rho  = r'*u; 
          beta = rho/rho_old; 
          q = u + beta*q; 
       end
       rho_old = rho; 
       tau_old = tau; 
       theta_old = theta; 
   end
   %%if (iter == maxit); solve_ok = 0; end; 
%%*************************************************************************
%% matvec: matrix-vector multiply.
%% matrix = [A.mat11 A.mat12; A.mat12' A.mat22]
%%*************************************************************************

   function Ax = matvec(A,x);

   m = length(A.mat11); m2 = length(x)-m; 
   if (m2 > 0)
      x1 = full(x(1:m)); 
   else
      x1 = full(x); 
   end
   Ax = mexMatvec(A.mat11,x1);
   if (m2 > 0)
      x2 = full(x(m+[1:m2]));
      Ax = Ax + mexMatvec(A.mat12,x2); 
      Ax2 = mexMatvec(A.mat12,x1,1) + mexMatvec(A.mat22,x2);
      Ax = [Ax; Ax2];  
   end
   return;
%%*************************************************************************
%% precond: 
%%*************************************************************************

   function Mx = precond(A,L,x)

   m = length(L.perm); m2 = length(x)-m;
   if (m2 > 0)
      x1 = full(x(1:m)); 
   else
      x1 = full(x); 
   end
   if (m2 > 0)
      x2 = x(m+[1:m2]);
      w = linsysolvefun(L,x1); 
      z = mexMatvec(A.mat12,w,1) -x2;
      z = L.Mu \ (L.Ml \ (L.Mp*z));
      x1 = x1 - mexMatvec(A.mat12,z); 
   end
%% 
   Mx = linsysolvefun(L,x1);  
%%
   if (m2 > 0)
      Mx = [Mx; z];
   end
%%*************************************************************************

