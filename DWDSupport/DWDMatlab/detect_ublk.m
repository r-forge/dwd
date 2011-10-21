%%*******************************************************************
%% detect_ublk: search for implied free variables in linear
%%              block. 
%% [blk2,At2,C2,ublkinfo] = detect_ublk(blk,At,C); 
%% 
%% Important: blk is assumed to have only 1 linear block.
%%
%% i1,i2: indices corresponding to splitting of unrestricted varaibles
%% i3   : remaining indices in the linear block
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%*******************************************************************

   function [blk2,At2,C2,ublkinfo] = detect_ublk(blk,At,C); 

   randn('seed',0);   
   blk2 = blk; At2 = At; C2 = C; 
   numblk = size(blk,1);
   ublkinfo = cell(size(blk,1),3);  

   for p = 1:numblk
      pblk = blk(p,:);
      m = size(At{p},2);        
      if strcmp(pblk{1},'l')
         r = randn(1,m);
         stime = cputime;
         Ap = At{p}'; Cp = C{p};
	 ApTr = (r*Ap)';
	 [sApTr,perm] = sort(abs(ApTr));
	 idx0 = find(abs(diff(sApTr)) < 1e-14);
         i1 = []; i2 = [];
	 if ~isempty(idx0)
            n = pblk{2}; 
            i1 = perm(idx0); i2 = perm(idx0+1);
	    Api1 = Ap(:,i1);
	    Api2 = Ap(:,i2);
	    Cpi1 = Cp(i1)';
	    Cpi2 = Cp(i2)';
	    idxzr = find(abs(Cpi1+Cpi2) < 1e-14 & sum(abs(Api1+Api2),1) < 1e-14);
            if ~isempty(idxzr)
               i1 = i1(idxzr');
               i2 = i2(idxzr');
               blk2{p,1} = 'u'; 
               blk2{p,2} = length(i1); 
               At2{p} = Ap(:,i1)'; 
               C2{p}  = Cp(i1); 
               fprintf(' %2.0d linear variables from unrestricted variable.\n',...
                         2*length(i1)); 
               i3 = setdiff([1:n],union(i1,i2));
               if ~isempty(i3)
                  blk2{numblk+1,1} = 'l'; 
                  blk2{numblk+1,2} = length(i3); 
                  At2{numblk+1,1} = Ap(:,i3)'; 
                  C2{numblk+1,1}  = Cp(i3); 
               end
	       ublkinfo{p,1} = i1; ublkinfo{p,2} = i2; ublkinfo{p,3} = i3;
            end
         end
      end
   end
%%*******************************************************************
