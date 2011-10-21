%%***********************************************************************
%% validate: validate data
%%
%% [At,C,dim,numblk,X0,Z0] = validate(blk,At,C,b,X0,y0,Z0);
%% 
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%***********************************************************************

   function  [At,C,dim,numblk,X0,Z0] = validate(blk,At,C,b,X0,y0,Z0);

   global spdensity

   if isempty(spdensity); spdensity = 0.5; end; 
%%
   if ~iscell(blk); 
      error('validate: blk must be a cell array'); end; 
   if size(blk,2)~=2
      error('validate: blk must be a cell array with 2 columns');       
   end 
   if ~iscell(At) | ~iscell(C); 
      error('validate: At, C must be cell arrays'); end;
   if (min(size(At))~=1 | min(size(C))~=1)
      error('validate: cell array At, C can only have 1 column or row'); end;
   if (size(At,2) > size(At,1));
      At = At'; end;
   if (size(C,2) > size(C,1));
      C = C'; end;
   if (nargin == 7)
      if ~iscell(X0) | ~iscell(Z0); 
         error('validate: X0, Z0 must be cell arrays'); 
      end
      if (min(size(X0))~=1 | min(size(Z0))~=1); 
         error('validate: cell array X, Z can only have 1 column or row'); 
      end
      if (size(X0,2) > size(X0,1)); X0 = X0'; end;
      if (size(Z0,2) > size(Z0,1)); Z0 = Z0'; end;
   end
%%
%% 
   m = length(b);  
   for p=1:size(blk,1)
      pblk = blk(p,:); 
      n = sum(pblk{2}); 
      if strcmp(pblk{1},'q') | strcmp(pblk{1},'l') | strcmp(pblk{1},'u'); 
         if (size(C{p},2) ~= 1); 
            error(['validate: ',num2str(p),'-th block of C must be column vectors']); 
         end;
         if (size(C{p},1) ~= n); 
            error(['validate: blk and C are not compatible']); 
         end; 
         if (size(At{p}) == [m n] & m~=n); 
            At{p} = At{p}'; end
         if ~all(size(At{p}) == [n,m]); 
            error('validate: blk and At not compatible'); end;
         if ~issparse(At{p}); 
            At{p} = sparse(At{p}); end; 
         if (nnz(C{p}) < spdensity*n); 
            if ~issparse(C{p}); C{p} = sparse(C{p}); end; 
         else
            if issparse(C{p}); C{p} = full(C{p}); end;
         end;
         if (nargin == 7)
            if ~all([size(X0{p},2) size(Z0{p},2)]==1); 
               error(['validate: ',num2str(p),'-th block of X0,Z0 must be column vectors']);
            end
            if ~all([size(X0{p},1) size(Z0{p},1)]==n); 
               error(['validate: blk, and X0,Z0, are not compatible']); 
            end              
            if (nnz(X0{p}) < spdensity*n); 
               if ~issparse(X0{p}); X0{p} = sparse(X0{p}); end; 
            else
               if issparse(X0{p}); X0{p} = full(X0{p}); end;
            end
            if (nnz(Z0{p}) < spdensity*n); 
               if ~issparse(Z0{p}); Z0{p} = sparse(Z0{p}); end; 
            else
               if issparse(Z0{p}); Z0{p} = full(Z0{p}); end;
            end
            if strcmp(pblk{1},'u') 
               Z0{p} = sparse(n,1); 
            end
         end
      else
         error(' blk: some fields are not specified correctly'); 
      end
   end
%%
%%-----------------------------------------
%% problem dimension
%%-----------------------------------------
%%
   dim = zeros(1,3);  numblk = 0;  
   for p = 1:size(blk,1)
      pblk = blk(p,:);
      if strcmp(pblk{1},'q')
         dim(1) = dim(1) + sum(pblk{2}); 
         numblk = numblk + length(pblk{2});
         nn(p) = length(pblk{2}); 
      elseif strcmp(pblk{1},'l')
         dim(2) = dim(2) + sum(pblk{2}); 
         nn(p) = sum(pblk{2}); 
      elseif strcmp(pblk{1},'u')
         dim(3) = dim(3) + sum(pblk{2}); 
         nn(p) = sum(pblk{2}); 
      end
   end
%%
%%***********************************************************************
