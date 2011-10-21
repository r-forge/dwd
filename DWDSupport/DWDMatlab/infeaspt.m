%%********************************************************************
%% infeaspt: generate an initial point for sdp.m
%%
%%  [X0,y0,Z0] = infeaspt(blk,At,C,b,options,scalefac);
%%
%%  options = 1  if want X0,Z0 to be scaled identity matrices
%%          = 2  if want X0,Z0 to be scalefac*(identity matrices).
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%********************************************************************

   function [X0,y0,Z0] = infeaspt(blk,At,C,b,options,scalefac);
%%
   if (nargin < 5); options = 1; end;
   if (options == 1); scalefac = []; end;
   if (options == 2) & (nargin < 6); scalefac = 1000; end;
   if (scalefac <= 0); error('scalefac must a positive number'); end;
%%
   if ~iscell(At); At = {At}; end;
   if ~iscell(C);  C = {C}; end;
   X0 = cell(size(C)); Z0 = cell(size(C));
   m = length(b); 
%%
   [At,C] = validate(blk,At,C,b);
%%
   for p = 1:size(blk,1); 
      pblk = blk(p,:); 
      blktmp = pblk{2};
      n = length(C{p});
      y0 = zeros(m,1);
      b2 = 1 + abs(b');
      if (options == 1);
         if strcmp(pblk{1},'q');
            s = 1+[0 cumsum(blktmp)];
            len = length(blktmp);
            normC = 1+norm(C{p});
            normA = 1+sqrt(sum(At{p}.*At{p}));
            idenqX = zeros(sum(blktmp),1);
            idenqX(s(1:len)) = sqrt(blktmp');
            X0{p} = max(1,max(b2./normA)) *idenqX;
            idenqZ = zeros(sum(blktmp),1);
	        normax = max([normA,normC]); 
            idenqZ(s(1:len)) = max([sqrt(blktmp); normax*ones(1,len)])';
            Z0{p} = idenqZ;
         elseif strcmp(pblk{1},'l');
            normC = 1+norm(C{p});
            normA = 1+sqrt(sum(At{p}.*At{p}));
            X0{p} = max(1,max(b2./normA)) *ones(n,1);
            Z0{p} = max(1,max([normA,normC])/sqrt(n)) *ones(n,1);
         elseif strcmp(pblk{1},'u');
            X0{p} = sparse(n,1);
            Z0{p} = sparse(n,1);
         else
            error(' blk: some fields not specified correctly'); 
         end
      elseif (options == 2);
         if strcmp(pblk{1},'q');
            s = 1+[0 cumsum(blktmp)];
            len = length(blktmp);
            idenq = zeros(sum(blktmp),1);
            idenq(s(1:len)) = ones(len,1);
            X0{p} = scalefac*idenq;
            Z0{p} = scalefac*idenq;
         elseif strcmp(pblk{1},'l');
            X0{p} = scalefac*ones(n,1);
            Z0{p} = scalefac*ones(n,1);
         elseif strcmp(pblk{1},'u');
            X0{p} = sparse(n,1);
            Z0{p} = sparse(n,1);
         else
            error(' blk: some fields not specified correctly'); 
         end
      end
   end
%%********************************************************************



