function [W,beta] = sep(X,y,penalty,t);

% This subroutine is designed to calculate a linear discriminator
% into K classes based on the training data X,
% the classification vector y, and the violation cost penalty.
% It uses the MDWD method of H.H. Zhang and M.J. Todd, based on the
% DWD method of J.S. Marron and M.J. Todd,
% "Distance Weighted Discrimination," TR 1339,
% School of Operations Research and Industrial Engineering,
% Cornell University, Ithaca, New York, available at
% ftp://ftp.orie.cornell.edu/pub/techreps/TR1339.pdf.

% The user will need to get the SDPT3 optimization package, Version 3.02,
% available at http://www.math.nus.edu.sg/~mattohkc/sdpt3.html,
% and install it. We recommend that sep.m be run from the SDPT3
% directory.

% SDPT3 includes an m-file startup.m, which sets the path,
% default parameters, and global environment variables.
% If the user has his/her own startup routine, this may not
% perform the necessary tasks. Either the instructions
% in SDPT3's startup.m should be appended, or SDPT3's
% startup.m be renamed startupSDPT3.m, say, and this called
% at the appropriate place.

% The subroutine contains various commented-out statements, which
% can be reinstated by those users wishing to see the time taken in
% various steps, the numbers of variables and constraints, etc.

% Input:  X, a d x n matrix, whose columns contain the training
%            data;
%         n x 1 vector y, with y_j the class of vector x_j; and
%         scalar penalty, the cost for perturbing each residual.

% Output: W, a d x K matrix, and beta, a K x 1 vector,
%            such that x is classified in class k if the kth
%            component of the K-vector W'x + beta is largest.
%         resid, the values of (w_k'x + beta_k) - (w_j'x + beta_j)
%            for x a column of X with label k.
%         totalviolation, the total amount added to the residuals.
%         dualgap, the duality gap, a measure of the accuracy
%            of the solutions found.
%         flag, an indication of the success of the computation:
%             0, success;
%            -1, inaccurate solution;
%            -2, problem infeasible or unbounded.

%stime = cputime;

% Find the dimensions of the data.

[d,n] = size(X);
nn = size(y,1);
if (nn ~= n), error('The dimensions are incompatible.'), end;
K = max(y);
if (K < 2), error('Fewer than 2 classes.'), end;

% Do the dimension reduction if in HDLSS setting.

if (d > n),
    [Q,XX] = qr(X,0);
    %qrtime = cputime - stime, stime = cputime;
    dnew = n;
else,
    XX = X;
    dnew = d;
    Q = speye(d);
end;
% t = d;

% nresid is the number of residuals,
% nv is the number of variables, and
% nc the number of constraints.

nresid = n*(K-1);
nv = K + K*dnew + 1 + 3*nresid + nresid + 1 + 2*K*d;
nc = 1 + 2*nresid + dnew + 1 + 1 + K*d;
e = ones(nresid,1);
speyen = speye(nresid);
%nv,
%nc,
% Set up the block structure, constraint matrix, rhs, and cost vector.

blk = cell(3,2);
blk{1,1} = 'u';
blk{1,2} = K;
blk{2,1} = 'q';
blk{2,2} = [K*dnew+1, 3*e'];
blk{3,1} = 'l';
blk{3,2} = nresid + 1 + 2*K*d;

Avec = cell(3,1);
A1 = sparse(nc,K);
A2 = sparse(nc,1+K*dnew+3*nresid);
A3 = sparse(nc,nresid);
A4 = sparse(nc,1+2*K*d);
rind = 0;
rin2 = nresid+1;
rin3 = 1+2*nresid + 1 + dnew;
cind = 1 + K*dnew;
for k = 1:K,
    indk = find(y == k);
    if isempty(indk), error('Empty class.'), end;
    Xk = XX(:,indk);
    numk = length(indk);
    for j = 1:K,
        if (j ~= k),
            A1(rind+1:rind+numk,j) = -ones(numk,1);
            A1(rind+1:rind+numk,k) =  ones(numk,1);
            A2(rind+1:rind+numk,(j-1)*dnew+2:j*dnew+1) = -Xk';
            A2(rind+1:rind+numk,(k-1)*dnew+2:k*dnew+1) =  Xk';
            A2(rind+1:rind+numk,cind+1:3:cind+0+3*numk) = -speye(numk);
            A2(rind+1:rind+numk,cind+2:3:cind+1+3*numk) =  speye(numk);
            A2(rin2+1:rin2+numk,cind+3:3:cind+2+3*numk) =  speye(numk);
            rind = rind + numk;
            rin2 = rin2 + numk;
            cind = cind + 3*numk;
        end;
    end;
    A2(1+2*nresid+2:1+2*nresid+1+dnew,(k-1)*dnew+2:k*dnew+1) = speye(dnew);
    A2(rin3+(k-1)*d+2:rin3+k*d+1,(k-1)*dnew+2:k*dnew+1) = Q;
end;
A1(1+2*nresid+1,:) = ones(1,K);
A2(nresid+1,1) = 1;
A3(1:nresid,:) = speyen;
A4(rin3+1,:) = ones(2*K*d+1,1);
A4(rin3+2:end,2:K*d+1) = speye(K*d);
A4(rin3+2:end,K*d+2:2*K*d+1) = -speye(K*d);
Avec{1,1} = A1';
Avec{2,1} = A2';
Avec{3,1} = [A3,A4]';  %'
%full(A1),
%full(A2),

b = [zeros(nresid,1);1;e;0;zeros(dnew,1);t;zeros(K*d,1)];
%b,

weighted_index = 1;
C = cell(3,1);
meann = n/K;
if weighted_index==0;
    c = zeros(nv-nresid-K-1-2*K*d,1);
    c(K*dnew+2:3:K*dnew+1+3*nresid) = e;
    c(K*dnew+3:3:K*dnew+2+3*nresid) = e;
    %c,
    C{1,1} = zeros(K,1);
    C{2,1} = c;
    C{3,1} = [penalty*e;zeros(1+2*K*d,1)];
else
    weight = ones(K,1);
    c = ones(nv-nresid-K-1-2*K*d,1);
    c1 = ones(nresid,1);
    rind = K*dnew+1;
    cind = 0;
    for k = 1:K,
        indk = find(y == k);
        if isempty(indk), error('Empty class.'), end;
        numk = length(indk);
        weight(k) = meann/numk;
        c(rind+1:3:rind-2+3*(K-1)*numk) = weight(k)*ones((K-1)*numk,1);
        c(rind+2:3:rind-1+3*(K-1)*numk) = weight(k)*ones((K-1)*numk,1);
        c1(cind+1:cind+(K-1)*numk) = weight(k)*penalty*ones((K-1)*numk,1);
        rind = rind + 3*(K-1)*numk;
        cind = cind + (K-1)*numk;
    end;
    C{1,1} = zeros(K,1);
    C{2,1} = c;
    C{3,1} = [c1;zeros(1+2*K*d,1)];
end;
%setuptime = cputime - stime, stime = cputime;

% Solve the SOCP problem.

% startup;
%startupSDPT3;
% Due to modification of parameters, variable OPTIONS is needed -- EZ
% parameters ;
[OPTIONS] = sqlparameters ;
OPTIONS.maxit = 40;
OPTIONS.vers = 2;
OPTIONS.steptol = 1e-10;
%OPTIONS.scale_data = 1,
OPTIONS.gaptol = 1e-12;
[X0,lambda0,Z0] = infeaspt(blk,Avec,C,b);
[obj,X,lambda,Z,info] = sqlp(blk,Avec,C,b,OPTIONS,X0,lambda0,Z0);
% If infeasible or unbounded, break.

if (info.termcode > 0), flag = -2; return; end;

% Compute the normal vectors W and constant terms beta.
%pause;

beta = X{1}; X2 = X{2}; X3 = X{3};
barW = reshape(X2(2:K*dnew+1),dnew,K);
if (d>n),
    W = Q*barW;
else,
    W = barW;
end;
return
