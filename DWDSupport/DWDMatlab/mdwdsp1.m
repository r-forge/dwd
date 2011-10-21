function [dirmat,beta] = MDWD(X,y,DWDpar,t) ;
% MDWD, Multiclass Distance Weighted Discrimination DIRection Matrix
%   Hanwen Huang matlab function
%
% Inputs:
%     X      - a d x n matrix, whose columns contain the training data;
%     y      - n x 1 vector, with y_j the class of vector x_j
%     DWDpar - the scalar penalty,
%                   (will adjust by median pairwise distance,
%                    100 is default)
%
% Output:
%    dirmat - direction matrix
%

%    Copyright (c) Hanwen Huang 2010



if nargin > 2 ;    %  then have input a threshfact, so use it
    threshfact = DWDpar ;
else ;    %  then use default threshfact
    threshfact = 100 ;
end ;



global CACHE_SIZE   % cache size in kbytes
global LOOP_LEVEL   % loop unrolling level
CACHE_SIZE = 256;
LOOP_LEVEL = 8;
%% set global variables for functions imported from LIPSOL


%  Compute median of pairwise distances squared between classes
%
n = size(X,2);
K = length(unique(y));
vpwdist2 = [] ;
for j = 1:(K-1);
    nj = length(find(y==j));
    xj = X(:,find(y==j));
    for k = (j+1):K;
        nk = length(find(y==k));
        xk = X(:,find(y==k));
        for i = 1:nj;
            pwdist2 = sum((vec2matSM(xj(:,i),nk) - xk).^2,1);
            vpwdist2 = [vpwdist2 pwdist2] ;
        end;
    end;
end;

medianpwdist2 = median(vpwdist2) ;

penalty = threshfact/medianpwdist2 ;
disp(penalty);
%  threshfact "makes this large",
%  and 1/medianpwdist2 "puts on correct scale"
if nargin<4;
    Wm = sepmdwd(X,y,penalty);
    t = sum(sum(abs(Wm)));
end;
[W,beta] = sepmdwdsp1(X,y,penalty,t);


% dirmat = W/sqrt(sum(sum(W.^2)));
% beta = beta/sqrt(sum(sum(W.^2)));
dirmat = W;
beta = beta;
