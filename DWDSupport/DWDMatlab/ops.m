%%******************************************************************
%% ops:  
%%
%%   Z = ops(X,operand,Y,alpha); 
%%
%%  INPUT:        X = a matrix or a scalar
%%                    or a CELL ARRAY consisting only of matrices 
%%          operand = sym, transpose, triu, tril,
%%                    real, imag, sqrt, abs, max, min, nnz,
%%                    spdiags, ones, norm, sum, row-norm, 
%%                    rank1, rank1inv, inv
%%                    +,  -, *, .*,  ./, .^ 
%%     Y = empty or a matrix or a scalar 
%%               or a CELL ARRAY consisting only of matrices
%%    alpha  = empty or a scalar 
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%******************************************************************

   function Z = ops(X,operand,Y,alpha); 

   if (nargin==2)
      if strcmp(operand,'norm')
         Z = 0; 
         for p = 1:length(X); Z = Z + sum(sum(X{p}.*X{p})); end;
         Z = sqrt(Z); 
      elseif strcmp(operand,'getM')
         Z = 0;
         for p = 1:length(X); Z = Z + size(X{p},1); end;
      end
   elseif (nargin==3)
      if strcmp(operand,'+') | strcmp(operand,'-') | strcmp(operand,'.*') ; 
         Z = cell(size(X)); 
         if strcmp(operand,'+')
            for p = 1:length(X); Z{p} = X{p} + Y{p}; end
         elseif strcmp(operand,'-')
            for p = 1:length(X); Z{p} = X{p} - Y{p}; end
         elseif strcmp(operand,'.*')
            for p = 1:length(X); Z{p} = X{p} .* Y{p}; end
         end
      else
         error([operand,' is not available']); 
      end
   else
      Z = cell(size(X)); 
      if strcmp(operand,'+') 
         for p = 1:length(X); Z{p} = X{p} + alpha*Y{p}; end
      elseif strcmp(operand,'-') 
         for p = 1:length(X); Z{p} = X{p} - alpha*Y{p}; end
      else
         error([operand,' is not available']); 
      end
   end
%%******************************************************************


