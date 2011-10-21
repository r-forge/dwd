function [tstat,epval,gfpval,zscore] = DiProPermSM(mdata1,mdata2,paramstruct) ;
% DiProPermSM, DIrection PROjection PERMutation based
% mean hypothesis test, intended for High Dimension, 
% Low Sample Size settings
%   Steve Marron's matlab function
%     Idea is to provide a 2 sample Mean Hypothesis Test,
%     based on a specified "directions" (e.g. DWD) in
%     High Dimension Low Sample Size contexts.
%     Method starts with an "important direction 
%     separating the two populations", then projects both
%     populations onto the direction vector, and next 
%     constructs a two sample t-statistic.  
%     Statistical significance is assessed, using a 
%     permutation approach, where the combined data is 
%     randomly relabelled into classes (and the direction
%     is recomputed for each permutation).
%
% Inputs:
%     mdata1 - d x n1 data set 1
%     mdata2 - d x n2 data set 2
%         i.e. data vectors are columns, with same nuber of rows
%
%     paramstruct - a Matlab structure of input parameters
%                      Use: "help struct" and "help datatypes" to
%                           learn about these.
%                      Create one, using commands of the form:
%
%       paramstruct = struct('field1',values1,...
%                            'field2',values2,...
%                            'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%
%    fields            values
%
%    idir             index for direction vector
%                         1  (default) DWD direction vector
%                         2  Mean Difference (aka centroid) direction vector
%                         3  Maximal Data Piling direction vector
%                         4  Fisher Linear Discrimination direction vector
%                                (based on generalized inverse)
%                         5  Support Vector Machine direction vector
%
%    ipval            index for p-value computation:
%                         1  (default) Display both empirical quantile
%                                and Gaussian fit quantile
%                         2  Show only emprical quantile
%                         3  Show only Gaussian fit quantile
%                         4  Show empirical quantile, Gaussian fit 
%                                quantile and also Z-score corresponding 
%                                to Gaussian fit quantile (allows for
%                                comparisons when quantiles are all are 0)
%
%    nsim             Number of simulated relabellings to use
%                           (default is 1000)
%                     0  -  only do DWD projection, and computation of t-statistic
%                               with no permutation test for signficance.
%                               (only gives the first graphical output,
%                                returns epval,gfpval,zscore as empty)
%
%    nreport          How often to report step taken to screen
%                     default 100 (report after 100 steps)
%                     (only has effect when iscreenwrite = 1)
%
%    seed             Seed for random number generator, use to allow
%                     identical repetition of random relabellings.
%                     Empty (default) for using the current value
%
%    icolor           0  fully black and white version (everywhere)
%                     1  (default)  color version (Red for Class 1, Blue for Class 2)
%                     2x3 color matrix, top row for Class 1, bottom row for Class 2
%                           Note:  larger color matrices are deliberately not 
%                                  supported, since color provides important
%                                  cues as to what the two classes are.
%                                  To combine multiple classes (having the
%                                  same color), consider a formulation like:
%                                       unique(icolor,'rows')
%                                
%
%    markerstr        Can be any of:
%                         A single string with symbol to use for marker,
%                             e.g. 'o' (default), '.', '+', 'x'
%                             (see "help plot" for a full list)
%                         A character array (2 x 1), of these symbols,
%                             One for each data set, created using:  strvcat
%                                                     or using:  ['+';'o']
%                         An (n1 + n2) x 1 character array, 
%                             One for each data point
%
%    ibigdot          0  (default)  use Matlab default for dot sizes
%                     1  force large dot size in prints (useful since some
%                              postscript graphics leave dots too small)
%                              (Caution: shows up as small in Matlab view)
%                              Only has effect on plot of simulated t-stats 
%
%    datovlaymax      maximum (on [0,1] scale, with 0 at bottom, 1 at top of plot)
%                     of vertical range for overlaid data.  Default = 0.6
%                     This applies to both the projected data, and the simulated t-stats
%
%    datovlaymin      minimum (on [0,1] scale, with 0 at bottom, 1 at top of plot)
%                     of vertical range for overlaid data.  Default = 0.5
%                     This applies to both the projected data, and the simulated t-stats
%
%    legendcellstr    cell array of strings for data set labels (2 of them),
%                     useful for (colored) classes, create this using
%                     cellstr, or {{string1 string2}}
%                     Also can indicate symbols, by just adding (at least 
%                             for +,x.o) into the text
%
%    title1str         string with title for left hand plot (showing projected data)
%                           default is empty string, '', to use:
%                               'Projections on DWD Direction' 
%                           for no title, use ' '
%
%    title2str         string with title for right hand plot (showing simulations)
%                           default is empty string, '', to use:
%                               [num2str(nsim) ' t-stats, from random relab''s'] 
%                           for no title, use ' '
%
%    titlefontsize    font size for title
%                                    (only has effect when the titlestr is nonempty)
%                           default is empty [], for Matlab default
%
%    xlabel1str        string with x axis label for left hand plot
%                           default is empty string, '', for no xlabel
%
%    xlabel2str        string with x axis label for left hand plot
%                           default is empty string, '', for no xlabel
%
%    ylabel1str        string with y axis label for left hand plot
%                           default is empty string, '', for no ylabel
%
%    ylabel2str        string with y axis label for left hand plot
%                           default is empty string, '', for no ylabel
%
%    labelfontsize     font size for axis labels
%                                    (only has effect when plot is made here,
%                                     and when a label str is nonempty)
%                           default is empty [], for Matlab default
%
%    vaxh              vector of 2 axis handles, useful for putting these two plots
%                      into chosen subplots.  E.g. can do:
%                          axh1 = subplot(2,2,2) ;
%                          axh2 = subplot(2,2,4) ;
%                          vaxh = [axh1; axh2] ;
%
%    DWDpar          DWDpenalty factor
%                         (will adjust by median pairwise dist.
%                              100 is default)
%                         (inactive, unless idir = 1)
%
%    SVMpar          SVMpenalty factor
%                         (SVMpar = C, when >= 0,
%                          SVMpar = penalty factor, when < 0
%                              (will adjust by median pairwise dist.)
%                         SVMpar = 1000 is default
%                         (inactive, unless idir = 5)
%
%    savestr          string controlling saving of output,
%                         either a full path, or a file prefix to
%                         save in matlab's current directory
%                         Will add .ps, and save as either
%                             color postscript (icolor ~= 0)
%                         or
%                             black&white postscript (when icolor = 0)
%                         unspecified:  results only appear on screen
%
%    iscreenwrite     0  (default)  no screen writes
%                     1  write to screen to show progress
%     
% Outputs:
%     Graphics in current Figure, 
%         showing projected data in left hand plot
%         and population of simulated pvalues in right hand plot
%     When savestr exists,
%        Postscript files saved in 'savestr'.ps
%                 (color postscript for icolor ~= 0)
%                 (B & W postscript for icolor = 0)
%     
%      tstat - 2 sample t statistic, computed for data projected 
%                  onto DWD direction vector
%
%      epval - empirical pvalue, based on simulated quantiles.
%                 summarizing results of permutation test
%
%     gfpval - Gaussian fit pvalue, based on Gaussian fit quantiles.
%                 summarizing results of permutation test
%
%     zscore - Z-score summary of gfpval, useful for comparisons, when
%                 gfpval = 0
%


% Assumes path can find personal functions:
%    DWD1SM.m
%    SVM1SM.m
%    projplot1SM.m
%    cprobSM.m
%    axisSM.m
%    kdeSM.m
%    nmfSM.m
%    bwsjpiSM.m
%    lbinrSM.m
%    vec2matSM.m
%    bwrfphSM.m
%    bwosSM.m
%    rootfSM
%    vec2matSM.m
%    bwrotSM.m
%    bwsnrSM.m
%    iqrSM.m
%    cquantSM.m

%    Copyright (c) J. S. Marron 2005



%  First set all parameters to defaults
%
idir = 1 ;
ipval = 1 ;
nsim = 1000 ;
nreport = 100 ;
seed = [] ;
icolor = 1 ;
markerstr = 'o' ;
ibigdot = 0 ;
idatovlay = 1 ;
datovlaymax = 0.6 ;
datovlaymin = 0.5 ;
legendcellstr = {} ;
title1str = '' ;
title2str = '' ;
titlefontsize = [] ;
xlabel1str = '' ;
xlabel2str = '' ;
ylabel1str = '' ;
ylabel2str = '' ;
labelfontsize = [] ;
vaxh = [] ;
DWDpar = 100 ;
SVMpar = 1000 ;
savestr = [] ;
iscreenwrite = 0 ;



%  Now update parameters as specified,
%  by parameter structure (if it is used)
%
if nargin > 2 ;   %  then paramstruct is an argument

  if isfield(paramstruct,'idir') ;    %  then change to input value
    idir = getfield(paramstruct,'idir') ; 
  end ;

  if isfield(paramstruct,'ipval') ;    %  then change to input value
    ipval = getfield(paramstruct,'ipval') ; 
  end ;

  if isfield(paramstruct,'nsim') ;    %  then change to input value
    nsim = getfield(paramstruct,'nsim') ; 
  end ;

  if isfield(paramstruct,'nreport') ;    %  then change to input value
    nreport = getfield(paramstruct,'nreport') ; 
  end ;

  if isfield(paramstruct,'seed') ;    %  then change to input value
    seed = getfield(paramstruct,'seed') ; 
  end ;

  if isfield(paramstruct,'icolor') ;    %  then change to input value
    icolor = getfield(paramstruct,'icolor') ; 
  end ;

  if isfield(paramstruct,'markerstr') ;    %  then change to input value
    markerstr = getfield(paramstruct,'markerstr') ; 
  end ;

  if isfield(paramstruct,'ibigdot') ;    %  then change to input value
    ibigdot = getfield(paramstruct,'ibigdot') ; 
  end ;

  if isfield(paramstruct,'idatovlay') ;    %  then change to input value
    idatovlay = getfield(paramstruct,'idatovlay') ; 
  end ;

  if isfield(paramstruct,'datovlaymax') ;    %  then change to input value
    datovlaymax = getfield(paramstruct,'datovlaymax') ; 
  end ;

  if isfield(paramstruct,'datovlaymin') ;    %  then change to input value
    datovlaymin = getfield(paramstruct,'datovlaymin') ; 
  end ;

  if isfield(paramstruct,'legendcellstr') ;    %  then change to input value
    legendcellstr = getfield(paramstruct,'legendcellstr') ; 
  end ;

  if isfield(paramstruct,'title1str') ;    %  then change to input value
    title1str = getfield(paramstruct,'title1str') ; 
  end ;

  if isfield(paramstruct,'title2str') ;    %  then change to input value
    title2str = getfield(paramstruct,'title2str') ; 
  end ;

  if isfield(paramstruct,'titlefontsize') ;    %  then change to input value
    titlefontsize = getfield(paramstruct,'titlefontsize') ; 
  end ;

  if isfield(paramstruct,'xlabel1str') ;    %  then change to input value
    xlabel1str = getfield(paramstruct,'xlabel1str') ; 
  end ;

  if isfield(paramstruct,'xlabel2str') ;    %  then change to input value
    xlabel2str = getfield(paramstruct,'xlabel2str') ; 
  end ;

  if isfield(paramstruct,'ylabel1str') ;    %  then change to input value
    ylabel1str = getfield(paramstruct,'ylabel1str') ; 
  end ;

  if isfield(paramstruct,'ylabel2str') ;    %  then change to input value
    ylabel2str = getfield(paramstruct,'ylabel2str') ; 
  end ;

  if isfield(paramstruct,'labelfontsize') ;    %  then change to input value
    labelfontsize = getfield(paramstruct,'labelfontsize') ; 
  end ;

  if isfield(paramstruct,'vaxh') ;    %  then change to input value
    vaxh = getfield(paramstruct,'vaxh') ; 
    if ~isempty(vaxh) ;
      if ~(sum(ishandle(vaxh)) == 2) ;
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
        disp('!!!   Warning from DiProPermSM.m:     !!!') ;
        disp('!!!   Invalid vaxh,                   !!!') ;
        disp('!!!   using default of two plots      !!!') ;
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
        vaxh = [] ;
      end ;
    end ;
  end ;

  if isfield(paramstruct,'DWDpar') ;    %  then change to input value
    DWDpar = getfield(paramstruct,'DWDpar') ; 
  end ;

  if isfield(paramstruct,'SVMpar') ;    %  then change to input value
    SVMpar = getfield(paramstruct,'SVMpar') ; 
  end ;

  if isfield(paramstruct,'savestr') ;    %  then use input value
    savestr = getfield(paramstruct,'savestr') ; 
    if ~(ischar(savestr) | isempty(savestr)) ;    %  then invalid input, so give warning
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from DiProPermSM.m:    !!!') ;
      disp('!!!   Invalid savestr,               !!!') ;
      disp('!!!   using default of no save       !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      savestr = [] ;
    end ;
  end ;

  if isfield(paramstruct,'iscreenwrite') ;    %  then change to input value
    iscreenwrite = getfield(paramstruct,'iscreenwrite') ; 
  end ;


end ;    %  of resetting of input parameters



%  Initiate parameters
%
d = size(mdata1,1) ;
if ~(d == size(mdata2,1)) ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Error from DiProPermSM.m:        !!!') ;
  disp('!!!   mdata1 and mdata2 must have      !!!') ;
  disp('!!!   same number of rows              !!!') ;
  disp('!!!   Terminating Execution            !!!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  tstat = [] ;
  epval = [] ;
  gfpval = [] ;
  zscore = [] ;
  return ;
end ;
n1 = size(mdata1,2) ;
n2 = size(mdata2,2) ;
n = n1 + n2 ;
mdata = [mdata1 mdata2] ;

if ~isempty(seed) ;
  rand('state',seed) ;
end ;

if nsim == 0 ;
  nax = 1 ;
else ;
  nax = 2 ;
end ;

if icolor == 0 ;    %  fully black and white version (everywhere)
  mcolor = zeros(n1 + n2,3) ;
  mlegcol = zeros(2,3) ;
  statstrcol = 'k' ;
elseif icolor == 1 ;    %  (default)  color version (Red for Class 1, Blue for Class 2)
  mcolor = [ones(n1,1) * [1 0 0] ; ...
            ones(n2,1) * [0 0 1]] ;
  mlegcol = [[1 0 0]; [0 0 1]] ;
  statstrcol = 'g' ;
elseif  size(icolor,1) == 2  &  size(icolor,2) == 3  ;
                  %  have 2 rows two icolor, as required
  mcolor = [ones(n1,1) * icolor(1,:) ; ...
            ones(n2,1) * icolor(2,:)] ;
  mlegcol = icolor ;
  statstrcol = 'g' ;
else ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Warning from DiProPermSM:        !!!') ;
  disp('!!!   Input icolor is invalid          !!!') ;
  disp(['!!!   Need 2 rows,    input has ' num2str(size(icolor,1))]) ;
  disp(['!!!   Need 3 columns, input has ' num2str(size(icolor,2))]) ;
  disp('!!!   Replacing icolor with default    !!!') ;
  disp('!!!       and proceeding               !!!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  mcolor = [ones(n1,1) * [1 0 0] ; ...
            ones(n2,1) * [0 0 1]] ;
  mlegcol = [[1 0 0]; [0 0 1]] ;
  statstrcol = 'g' ;
end ;

if size(markerstr,1) == 1 ;    %  then can use as single symbol
  vmarkerstr = markerstr ;
elseif size(markerstr,1) == 2 ;    %  then need to expand out
  vmarkerstr = [] ;
  for i = 1:n1 ;
    vmarkerstr = strvcat(vmarkerstr,markerstr(1,1)) ;
  end ;
  for i = 1:n2 ;
    vmarkerstr = strvcat(vmarkerstr,markerstr(2,1)) ;
  end ;
else ;    %  pass given markerstr through to graphics
  vmarkerstr = markerstr ;
end ;



%  Find Data 1 vs. Data 2 direction vector
%
if idir == 1 ;    %  DWD direction vector

  dirstr = 'DWD' ;
  vdir = DWD1SM(mdata1,mdata2,DWDpar) ;
      %  DWD direction vector, pointing from 2nd group towards first

elseif idir == 2 ;    %  Mean Difference (aka centroid) direction vector

  dirstr = 'MD' ;
  vdir = mean(mdata1,2) - mean(mdata2,2) ;
      %  MD direction vector, pointing from 2nd group towards first

elseif idir == 3 ;    %  Maximal Data Piling direction vector

  dirstr = 'MDP' ;
  vmean1 = mean(mdata1,2) ;
  vmean2 = mean(mdata2,2) ;
  mtrain = [mdata1 mdata2] ;
  vcent = mean(mtrain,2) ;
      %  centerpoint
  mresid = mtrain - vec2matSM(vcent,size(mtrain,2)) ;
  mcov = cov(mresid') ;
      %  Get covariance matrix, transpose, since want 
      %               "coordinates as variables"
      %  This gives "global covariance"
  mcovinv = pinv(mcov) ;
      %  pseudo-inverse
  vdir = mcovinv * (vmean1 - vmean2) ;
      %  Maximal Data Piling Vector

elseif idir == 4 ;    %  Fisher Linear Discrimination direction vector

  dirstr = 'FLD' ;
  vmean1 = mean(mdata1,2) ;
  vmean2 = mean(mdata2,2) ;
  mresid1 = mdata1 - vec2matSM(vmean1,size(mdata1,2)) ;
  mresid2 = mdata2 - vec2matSM(vmean2,size(mdata2,2)) ;
  mresid = [mresid1 mresid2] ;
  mcov = cov(mresid') ;
      %  Get covariance matrix, transpose, since want 
      %               "coordinates as variables"
      %  This gives "pooled within class covariance"
  mcovinv = pinv(mcov) ;
      %  pseudo-inverse
  vdir = mcovinv * (vmean1 - vmean2) ;
      %  Fisher Linear Discriminant Vector

elseif idir == 5 ;    %  Support Vector Machine direction vector

  dirstr = 'SVM' ;
  vdir = SVM1SM(mdata1,mdata2,SVMpar) ;
      %  SVM direction vector, pointing from 2nd group towards first

else ;

  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Error from DiProPermSM.m:     !!!') ;
  disp('!!!       idir is invalid           !!!') ;
  disp('!!!   Terminating exceution         !!!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  tstat = [] ;
  epval = [] ;
  gfpval = [] ;
  zscore = [] ;
      %  set to empty, to avoid error.
  return ;

end ;
vdir = vdir / sqrt(sum(vdir.^2)) ;
    %  divide by length, to make it a unit vector



%  Make Hypo test output plot
%
if isempty(vaxh) ;
  subplot(1,nax,1) ;
else ;
  axes(vaxh(1)) ;
end ;

if isempty(title1str) ;
  projtitstr = ['Projections on ' dirstr ' Direction']  ;
else ;
  projtitstr = title1str ;
end ;

paramstructPP1 = struct('icolor',mcolor, ...
                        'isubpopkde',1, ...
                        'markerstr',vmarkerstr, ...
                        'titlestr',projtitstr, ...
                        'titlefontsize',titlefontsize, ...
                        'xlabelstr',xlabel1str, ...
                        'ylabelstr',ylabel1str, ...
                        'labelfontsize',labelfontsize, ...
                        'datovlaymin',datovlaymin, ...
                        'datovlaymax',datovlaymax, ...
                        'iscreenwrite',iscreenwrite) ;
if ~(isempty(legendcellstr)) ;
  paramstructPP1 = setfield(paramstructPP1,'legendcellstr',legendcellstr) ;
  paramstructPP1 = setfield(paramstructPP1,'mlegendcolor',mlegcol) ;
end ;

projplot1SM(mdata,vdir,paramstructPP1) ;

%  overlay combined hypo test t-stat
%
if idir ~= 3 ;    %  then are doing somethings besides MDP
  [h,pval,ci,stats] = ttest2(mdata1' * vdir, ...
                             mdata2' * vdir) ;
  tstat = getfield(stats,'tstat') ;
  tstatstr = 't-stat' ;
else ;    %  then are doing MPD, so don't compute T-stat
  tstat = mean(mdata1' * vdir) - mean(mdata2' * vdir) ;
  tstatstr = 'M-Diff' ;
end ;
vax = axis ;
hold on ;
  text(vax(1) + 0.1 * (vax(2) - vax(1)), ...
       vax(3) + 0.92 * (vax(4) - vax(3)), ...
       [tstatstr ' = ' num2str(tstat)], ...
       'Color',statstrcol) ;
hold off ;



if ~(nsim == 0) ;    %  Then do permutation test

  %  Recompute t-stats over random relabellings
  %
  vtstat = [] ;
  for isim = 1:nsim ;

    if (isim / nreport) == floor(isim / nreport) ;
      if iscreenwrite == 1 ;
        disp(['    Working on sim ' num2str(isim) ' of ' num2str(nsim)]) ;
      end ;
    end ;

    flagss1sim = [ones(1,n1), zeros(1,n2)] ;
    vunif = rand(1,n) ;
    [temp,indperm] = sort(vunif) ;
        %  indices of random permutation
    flagss1sim = flagss1sim(indperm) ;
        %  random permutation of flagss1sim

    comboss1flag = logical(flagss1sim) ; 


    if idir == 1 ;    %  DWD direction vector

      vdirc = DWD1SM(mdata(:,comboss1flag), ...
                       mdata(:,~comboss1flag)) ;
          %  DWD direction vector, pointing from 2nd group towards first

    elseif idir == 2 ;    %  Mean Difference (aka centroid) direction vector

      vdirc = mean(mdata(:,comboss1flag),2) - ...
                  mean(mdata(:,~comboss1flag),2) ;
          %  MD direction vector, pointing from 2nd group towards first

    elseif idir == 3 ;    %  Maximal Data Piling direction vector

      vmean1c = mean(mdata(:,comboss1flag),2) ;
      vmean2c = mean(mdata(:,~comboss1flag),2) ;
      mtrainc = [mdata(:,comboss1flag) mdata(:,~comboss1flag)] ;
      vcentc = mean(mtrainc,2) ;
          %  centerpoint
      mresidc = mtrainc - vec2matSM(vcentc,size(mtrainc,2)) ;
      mcovc = cov(mresidc') ;
          %  Get covariance matrix, transpose, since want 
          %               "coordinates as variables"
          %  This gives "global covariance"
      mcovinvc = pinv(mcovc) ;
          %  pseudo-inverse
      vdirc = mcovinvc * (vmean1c - vmean2c) ;
          %  Maximal Data Piling Vector

    elseif idir == 4 ;    %  Fisher Linear Discrimination direction vector

      vmean1c = mean(mdata(:,comboss1flag),2) ;
      vmean2c = mean(mdata(:,~comboss1flag),2) ;
      mresid1c = mdata(:,comboss1flag) - vec2matSM(vmean1c,sum(comboss1flag)) ;
      mresid2c = mdata(:,~comboss1flag) - vec2matSM(vmean2c,sum(~comboss1flag)) ;
      mresidc = [mresid1c mresid2c] ;
      mcovc = cov(mresidc') ;
          %  Get covariance matrix, transpose, since want 
          %               "coordinates as variables"
          %  This gives "pooled within class covariance"
      mcovinvc = pinv(mcovc) ;
          %  pseudo-inverse
      vdirc = mcovinvc * (vmean1c - vmean2c) ;
          %  Fisher Linear Discriminant Vector

    elseif idir == 5 ;    %  Support Vector Machine direction vector

      vdirc = SVM1SM(mdata(:,comboss1flag), ...
                       mdata(:,~comboss1flag)) ;
          %  SVM direction vector, pointing from 2nd group towards first

    end ;
    vdirc = vdirc / sqrt(sum(vdirc.^2)) ;
        %  divide by length, to make it a unit vector


    if idir ~= 3 ;    %  then are doing something besides MDP
      [h,pval,ci,stats] = ttest2(mdata(:,comboss1flag)' * vdirc, ...
                                 mdata(:,~comboss1flag)' * vdirc) ;
      tstatsim = getfield(stats,'tstat') ;
    else ;    %  then are doing MPD, so don't compute T-stat
      tstatsim = mean(mdata(:,comboss1flag)' * vdirc) - ...
              mean(mdata(:,~comboss1flag)' * vdirc) ;
    end ;

    vtstat = [vtstat; tstatsim] ;


  end ;


  %  Compute p-values
  %
  epval = 1 - cprobSM(vtstat,tstat) ;
      %  emprical p-value
  simmean = mean(vtstat) ;
  simsd = std(vtstat) ;
  zscore = (tstat - simmean) / simsd ;
  gfpval = 1 - normcdf(zscore) ;
      %  Gaussian fit p-value


  %  add results to graphics
  %
  if isempty(vaxh) ;
    subplot(1,2,2) ;
  else ;
    axes(vaxh(2)) ;
  end ;

  vax = axisSM([vtstat; tstat]) ;
  if isempty(title2str) ;
    kdetitstr = [num2str(nsim) ' ' tstatstr 's, from random relab''s'] ;
  else ;
    kdetitstr = title2str ;
  end ;

  kdeparamstruct = struct('vxgrid',vax, ...
                          'linecolor','k', ...
                          'dolcolor','k', ...
                          'ibigdot',ibigdot, ...
                          'titlestr',kdetitstr, ...
                          'titlefontsize',titlefontsize, ...
                          'xlabelstr',xlabel2str, ...
                          'ylabelstr',ylabel2str, ...
                          'labelfontsize',labelfontsize, ...
                          'datovlaymin',datovlaymin, ...
                          'datovlaymax',datovlaymax, ...
                          'iscreenwrite',1) ;
  kdeSM(vtstat,kdeparamstruct) ;
  vax = axis ;
  hold on ;
    plot([tstat; tstat],[vax(3); vax(4)],'Color',statstrcol) ;
    text(vax(1) + 0.3 * (vax(2) - vax(1)), ...
         vax(3) + 0.9 * (vax(4) - vax(3)), ...
         [tstatstr ' = ' num2str(tstat)],'Color',statstrcol) ;
    if ipval == 1 ;
      text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
           vax(3) + 0.8 * (vax(4) - vax(3)), ...
           ['Emprical pval = ' num2str(epval)],'Color','k') ;
      text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
           vax(3) + 0.7 * (vax(4) - vax(3)), ...
           ['Gaussian fit pval = ' num2str(gfpval)],'Color','k') ;
   elseif ipval == 2 ;
      text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
           vax(3) + 0.8 * (vax(4) - vax(3)), ...
           ['Emprical pval = ' num2str(epval)],'Color','k') ;
    elseif ipval == 3 ;
      text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
           vax(3) + 0.8 * (vax(4) - vax(3)), ...
           ['Gaussian fit pval = ' num2str(gfpval)],'Color','k') ;
    elseif ipval == 4 ;
      text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
           vax(3) + 0.8 * (vax(4) - vax(3)), ...
           ['Emprical pval = ' num2str(epval)],'Color','k') ;
      text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
           vax(3) + 0.7 * (vax(4) - vax(3)), ...
           ['Gaussian fit pval = ' num2str(gfpval)],'Color','k') ;
      text(vax(1) + 0.2 * (vax(2) - vax(1)), ...
           vax(3) + 0.6 * (vax(4) - vax(3)), ...
           ['Gaussian fit Z-score = ' num2str(zscore)],'Color','k') ;
    else ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from DiProPermSM.m:         !!!') ;
      disp('!!!   Invalid ipval,                      !!!') ;
      disp('!!!   using default of both quantiles     !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      text(vax(1) + 0.3 * (vax(2) - vax(1)), ...
           vax(3) + 0.8 * (vax(4) - vax(3)), ...
           ['Emprical pval = ' num2str(epval)],'Color','k') ;
      text(vax(1) + 0.3 * (vax(2) - vax(1)), ...
           vax(3) + 0.7 * (vax(4) - vax(3)), ...
           ['Gaussian fit pval = ' num2str(gfpval)],'Color','k') ;
    end ;

    if  ipval == 1  |  ipval == 3  |  ipval == 4  ;    %  Add Gaussian density
      xgrid = linspace(vax(1),vax(2),401)' ;
      plot(xgrid,nmfSM(xgrid,simmean,simsd^2,1),'k') ;
    end ;

  hold off ;


else ;    %  for nsim = 0

  epval = [] ;
  gfpval = [] ;
  zscore = [] ;
      %  set to empty, to avoid error.

end ;



%  Save graphical output (if needed)
%
if ~isempty(savestr) ;   %  then create postscript file

  orient landscape ;

  if  (size(icolor,2) > 1.5)  |  (icolor ~= 0)  ;     %  then make color postscript
    print('-dpsc',savestr) ;
  else ;                %  then make black and white
    print('-dps',savestr) ;
  end ;

end ;




