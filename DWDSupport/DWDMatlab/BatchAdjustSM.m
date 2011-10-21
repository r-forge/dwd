function adjdata = BatchAdjustSM(rawdata,batchlabels,paramstruct) 
% BATCHADJUSTSM, of adjustment of "centerpoint" of subpopulations,
%   Steve Marron's matlab function
%     Intended for simple adjustment of subpopulation mean effects
%     in microarray data analysis.
%     Allows output of graphics, including "before" and "after"
%     PCA 2-d Draftsman's plots.
%     Algorithm uses DWD to find an effectrive direction for adjustment,
%     projects the data in this direction, and subtracts the 
%     subpopulation means.
%     This works pairwise, i.e. for a pair of subpopulations
%     For more subpopulations, apply this several times, to groups
%     of subpopulations.  Good choice of groups can be done by first
%     looking at the population structure, e.g. using a class colored
%     version of the PCA 2-d Draftsman's plot given by curvdatSM.
% Inputs:
%   rawdata     - d x n matrix of log gene expression data, 
%                          columns are cases, rows are genes 
%   batchlabels - 1 x n vector of vector of batch labels,
%                          for each case, must be +-1
%   paramstruct - a Matlab structure of input parameters
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%       paramstruct = struct('field1',values1, ...
%                            'field2',values2, ...
%                            'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%
%    fields            values
%
%    viplot           vector of zeros and ones, indicating which plots to make
%                         1st entry:  1 (default) makes "before"
%                                           PCA 2-d Draftsman's plot
%                         2nd entry:  1 (default) makes "projection plot"
%                                           showing DWD performance
%                         3rd entry:  1 makes "afterwards projection plot",
%                                           of DWD applied to adjusted data
%                                           (default is 0, no plot) 
%                         4th entry:  1 (default) makes "after"
%                                           PCA 2-d Draftsman's plot
%                             (use zeros(4,1) for no plots)
%
%    imeantype        0   (default) move both projected populations to 0
%                                  (sensible for cDNA and other 
%                                   differentially expressed data)
%                     +1  move -1 data set so projected mean is same as 
%                                    +1 data set
%                     -1  move +1 data set so projected mean is same as 
%                                    -1 data set
%
%
%    savestr          string controlling saving of output,
%                         either a full path, or a file prefix to
%                         save in matlab's current directory
%                     unspecified:  results only appear on screen
%                     result:  add various plot names (depending 
%                                 on viplot) and add .ps
%
%    titlestr         String for Title of Projection plots 
%                           (will add to this depending on plot)
%                           (leave empty for no title at all)
%                           (default is "Batch Adjustment")
%
%    titlefontsize    Font Size for titles (uses Matlab default)
%                                   (18 is "fairly large")
%
%    legcellstr       Legend Cell String
%                     Use to apply labels to batches
%                     E.g.    legcellstr = {{'Batch 1' 'Batch 2'}} ;
%                         (this "cell within a cell" structure seems
%                          needed to pass in a cell array of strings,
%                          for Matlab 6.0 & 7.0, may be different for 
%                          other Matlab versions)
%
%    iscreenwrite     0  (default)  no screen writes
%                     1  write to screen to show progress
%
%    minproj          left end of range in projection plots
%                               (default is output of axisSM)
%
%    maxproj          right end of range in projection plots
%                               (default is output of axisSM)
%
%    npc              Number of Principal Components to use in 
%                     draftsman's plot of 2-d PC projections
%                     (default = 6)
%
%
% Output:
%   adjdata     - d x n matrix of adjusted data, 
%                          columns are cases, rows are genes 
%    

%    Copyright (c) J. S. Marron 2003-2005




%  Set all input parameters to defaults
%
viplot = [1; 1; 0; 1] ;
imeantype = 0 ;
savestr = [] ;
titlestr = ['Batch Adjustment'] ;
titlefontsize = [] ;
legcellstr = {} ;
iscreenwrite = 0 ;
minproj = [] ;
maxproj = [] ;
npc = 6 ;



%  Now update parameters as specified,
%  by parameter structure (if it is used)
%
if nargin > 2 ;   %  then paramstruct is an argument

  if isfield(paramstruct,'viplot') ;    %  then change to input value
    viplot = getfield(paramstruct,'viplot') ; 
  end ;

  if isfield(paramstruct,'imeantype') ;    %  then change to input value
    imeantype = getfield(paramstruct,'imeantype') ; 
  end ;

  if isfield(paramstruct,'savestr') ;    %  then use input value
    savestr = getfield(paramstruct,'savestr') ; 
    if ~ischar(savestr) ;    %  then invalid input, so give warning
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from BatchAdjustSM.m:    !!!') ;
      disp('!!!   Invalid savestr,                 !!!') ;
      disp('!!!   using default of no save         !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      savestr = [] ;
    end ;
  end ;

  if isfield(paramstruct,'titlestr') ;    %  then change to input value
    titlestr = getfield(paramstruct,'titlestr') ; 
  end ;

  if isfield(paramstruct,'titlefontsize') ;    %  then change to input value
    titlefontsize = getfield(paramstruct,'titlefontsize') ; 
  end ;

  if isfield(paramstruct,'legcellstr') ;    %  then change to input value
    legcellstr = getfield(paramstruct,'legcellstr') ; 
  end ;

  if isfield(paramstruct,'iscreenwrite') ;    %  then change to input value
    iscreenwrite = getfield(paramstruct,'iscreenwrite') ; 
  end ;

  if isfield(paramstruct,'minproj') ;    %  then change to input value
    minproj = getfield(paramstruct,'minproj') ; 
  end ;

  if isfield(paramstruct,'maxproj') ;    %  then change to input value
    maxproj = getfield(paramstruct,'maxproj') ; 
  end ;

  if isfield(paramstruct,'npc') ;    %  then change to input value
    npc = getfield(paramstruct,'npc') ; 
  end ;

end ;    %  of resetting of input parameters



%  Set internal parameters
%
npixshift = 20 ;
    %  number of pixels to shift new figures by



%  Check inputs
%
d = size(rawdata,1) ;
n = size(rawdata,2) ;
errflag = logical(0) ;
if  (n ~= size(batchlabels,2))  | ...
    (1 ~= size(batchlabels,1))  ;

  errflag = logical(1) ;

  errstr = ['input "batchlabels" must be a row vector ' ...
                     'of length ' num2str(n)] ;

elseif  (sum(batchlabels == 1) + sum(batchlabels == -1))  ~=  n  ;

  errflag = logical(1) ;

  errstr = 'entries in "batchlabels" must all be +- 1' ;
  
end ;

if ~iscell(legcellstr) ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Warning from BatchAdjustSM:          !!!') ;
  disp('!!!   legcellstr is not a cell string      !!!') ;
  disp('!!!   Will turn off the Legend Cell String !!!') ;
  disp('!!!   i.e. the text with batch labels      !!!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  legcellstr = {} ;

else ;

  if length(legcellstr) > 2 ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from BatchAdjustSM:       !!!') ;
    disp('!!!   legcellstr is too big,            !!!') ;
    disp('!!!   Will use only first two entries   !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    legcellstr = legcellstr([1,2]) ;

  end ;

end ;



if errflag ;    %  then had a fatal input error

  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Error from BatchAdjustSM.m:                   !!!') ;
  disp(['!!!   ' errstr]) ;
  disp('!!!   Terminating execution, with an empty return   !!!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;

  adjdata = [] ;


else ;    %  inputs OK, so do serious work


  if length(viplot) ~= 4 ;

    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from BatchAdjustSM.m:   !!!') ;
    disp('!!!   Invalid size of viplot,         !!!') ;
    disp('!!!   reverting to default            !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;

    viplot = [1; 1; 0; 1] ;
  
  end ;

  viplot = logical(viplot) ;
      %  make sure this is logical, 
      %  for use in if statements below



  %  Do DWD Batch adjustment
  %
  flagp = (batchlabels == 1) ;
  flagn = (batchlabels == -1) ;


  %  Find DWD direction
  %
  dirvec = DWD1SM(rawdata(:,flagp),rawdata(:,flagn)) ;


  %  Project data
  %
  vprojp = rawdata(:,flagp)' * dirvec ;
  vprojn = rawdata(:,flagn)' * dirvec ;

  meanprojp = mean(vprojp) ;
  meanprojn = mean(vprojn) ;


  %  Do shift along direction vector
  %
  if imeantype == -1 ;    %  move +1 data, so projected mean is same as -1 data
    adjdata(:,flagp) = rawdata(:,flagp) - vec2matSM(meanprojp * dirvec,length(vprojp)) ...
                                        + vec2matSM(meanprojn * dirvec,length(vprojp)) ;
    adjdata(:,flagn) = rawdata(:,flagn) ;
  elseif imeantype == 1 ;    %  move -1 data, so projected mean is same as +1 data
    adjdata(:,flagp) = rawdata(:,flagp) ;
    adjdata(:,flagn) = rawdata(:,flagn) - vec2matSM(meanprojn * dirvec,length(vprojn)) ...
                                        + vec2matSM(meanprojp * dirvec,length(vprojn)) ;
  else ;    %  move both projected populations to 0
    adjdata(:,flagp) = rawdata(:,flagp) - vec2matSM(meanprojp * dirvec,length(vprojp)) ;
    adjdata(:,flagn) = rawdata(:,flagn) - vec2matSM(meanprojn * dirvec,length(vprojn)) ;
  end ;



  %  Make graphics (if needed)
  %
  fignum = 1 ;
  if viplot(1) ;    %  then make "before" 2-d scatterplot matrix

    figure(fignum) ;
    clf ;
    fignum = fignum + 1 ;
    set(gcf,'Name','PCA 2-d Projections, Before') ;

    mlegendcolorBatch = [[1 0 0] ; ...
                        [0 0 1]] ;

    mcolorBatch = [(batchlabels == 1)' zeros(n,2)]+...
                  [ zeros(n,2) (batchlabels == -1)'];


    paramstruct = struct('viout',[3], ...
                         'vipcplot',1:npc, ...
                         'icolor',mcolorBatch, ...
                         'iscreenwrite',iscreenwrite) ;

    if ~isempty(savestr) ;
      paramstruct = setfield(paramstruct, ...
                            'savestr',[savestr 'PCA2dBefore']) ;
    end ;

    if ~isempty(legcellstr) ;
      paramstruct = setfield(paramstruct, ...
                            'legendcellstr',legcellstr) ;
      paramstruct = setfield(paramstruct, ...
                            'mlegendcolor',mlegendcolorBatch) ;
    end ;

    curvdatSM(rawdata,paramstruct) ;


  end ;



  if viplot(2) ;    %  then make DWD projection plot

    if fignum == 1 ;    %  then this figure is first
      figure(1) ;
      clf ;
    else ;    %  then have an earlier figure, so shift this one
      vpos = get(gcf,'Position') ;
      vpos = [(vpos(1) + npixshift) ...
              (vpos(2) - npixshift) ...
               vpos(3) vpos(4)] ;
      figure(fignum) ;
      clf ;
      set(gcf,'Position',vpos) ;
    end ;
    set(gcf,'Name','DWD Projections') ;
    fignum = fignum + 1 ;

    %  compute kernel density estimates
    %
    h = bwsjpiSM([vprojp; vprojn]) ;
    vax = axisSM([vprojp; vprojn]) ;
    if ~isempty(minproj) ;
      vax(1) = minproj ;
    end ;
    if ~isempty(maxproj) ;
      vax(2) = maxproj ;
    end ;
      paramstruct = struct('vh',h, ...
                           'vxgrid',[vax(1) vax(2)]) ;
    [vkdep1, xgrid1] = kdeSM(vprojp,paramstruct) ;
    vkden1 = kdeSM(vprojn,paramstruct) ;


    %  Make Projection Plots
    %
    clf ;
    plot(xgrid1,vkdep1,'r-','LineWidth',2) ;
      axisSM(xgrid1,[vkdep1;vkden1]) ;
      vax = axis ;
      if ~isempty(titlestr) ;
        titstr = [titlestr ', DWD Projection'] ;
        if isempty(titlefontsize) ;
          title(titstr) ;
        else ;
          title(titstr,'FontSize',titlefontsize) ;
        end ;
      end ;

      hold on ;
        plot(xgrid1,vkden1,'b-','LineWidth',2) ;

        plotbottom = vax(3) ;
        plottop = vax(4) ;
        yrand = plotbottom + (0.7 + 0.1 * rand(length(vprojp),1)) ...
                                              * (plottop - plotbottom) ;
               %  y coords for jittering
        plot(vprojp,yrand,'r+') ;

        yrand = plotbottom + (0.5 + 0.1 * rand(length(vprojn),1)) ...
                                              * (plottop - plotbottom) ;
               %  y coords for jittering
        plot(vprojn,yrand,'bo') ;

        if ~isempty(legcellstr) ;
            tx = vax(1) + 0.1 * (vax(2) - vax(1)) ;
            ty = vax(3) + 0.9 * (vax(4) - vax(3)) ;
            voutstr = legcellstr{1} ;
          text(tx,ty,voutstr,'Color','r','FontSize',18) ;
            ty = vax(3) + 0.4 * (vax(4) - vax(3)) ;
            vnotoutstr = legcellstr{2} ;
          text(tx,ty,vnotoutstr,'Color','b','FontSize',18) ;
        end ;

      hold off ;


    %  Save as .ps
    %
    if ~isempty(savestr) ;

      if iscreenwrite == 1 ;
        disp('    BatchAdjustSM saving results') ;
      end ;

      orient landscape ;
      print('-dpsc2',[savestr 'DWDproj']) ;

      if iscreenwrite == 1 ;
        disp('    BatchAdjustSM finished save') ;
        disp('  ') ;
      end ;

    end ;


  end ;



  if viplot(3) ;    %  then repeat DWD, for checking, and make projection plot

    if fignum == 1 ;    %  then this figure is first
      figure(1) ;
      clf ;
    else ;    %  then have an earlier figure, so shift this one
      vpos = get(gcf,'Position') ;
      vpos = [(vpos(1) + npixshift) ...
              (vpos(2) - npixshift) ...
               vpos(3) vpos(4)] ;
      figure(fignum) ;
      clf ;
      set(gcf,'Position',vpos) ;
    end ;
    set(gcf,'Name','Check DWD Proj. for Adjusted Data') ;
    fignum = fignum + 1 ;

    %  Find DWD direction for adjusted data
    %
    dirvec = DWD1SM(adjdata(:,flagp),adjdata(:,flagn)) ;

    %  Project data
    %
    vreprojp = adjdata(:,flagp)' * dirvec ;
    vreprojn = adjdata(:,flagn)' * dirvec ;

    %  compute kernel density estimates
    %
    h = bwsjpiSM([vreprojp; vreprojn]) ;
    vax = axisSM([vreprojp; vreprojn]) ;
    if ~isempty(minproj) ;
      vax(1) = minproj ;
    end ;
    if ~isempty(maxproj) ;
      vax(2) = maxproj ;
    end ;
      paramstruct = struct('vh',h, ...
                           'vxgrid',[vax(1) vax(2)]) ;
    [vkdep2, xgrid2] = kdeSM(vreprojp,paramstruct) ;
    vkden2 = kdeSM(vreprojn,paramstruct) ;

    %  Make Projection Plots
    %
    clf ;
    plot(xgrid2,vkdep2,'r-','LineWidth',2) ;
      axisSM(xgrid2,[vkdep2;vkden2]) ;
      vax = axis ;
      if ~isempty(titlestr) ;
        titstr = [titlestr ', Check DWD Proj. for Adj''d Data'] ;
        if isempty(titlefontsize) ;
          title(titstr) ;
        else ;
          title(titstr,'FontSize',titlefontsize) ;
        end ;
      end ;

      hold on ;
        plot(xgrid2,vkden2,'b-','LineWidth',2) ;

        plotbottom = vax(3) ;
        plottop = vax(4) ;
        yrand = plotbottom + (0.7 + 0.1 * rand(length(vreprojp),1)) ...
                                              * (plottop - plotbottom) ;
               %  y coords for jittering
        plot(vreprojp,yrand,'r+') ;

        yrand = plotbottom + (0.5 + 0.1 * rand(length(vreprojn),1)) ...
                                              * (plottop - plotbottom) ;
               %  y coords for jittering
        plot(vreprojn,yrand,'bo') ;

        if ~isempty(legcellstr) ;
            tx = vax(1) + 0.1 * (vax(2) - vax(1)) ;
            ty = vax(3) + 0.9 * (vax(4) - vax(3)) ;
            voutstr = legcellstr{1} ;
          text(tx,ty,voutstr,'Color','r','FontSize',18) ;
            ty = vax(3) + 0.4 * (vax(4) - vax(3)) ;
            vnotoutstr = legcellstr{2} ;
          text(tx,ty,vnotoutstr,'Color','b','FontSize',18) ;
        end ;

      hold off ;


    %  Save as .ps
    %
    if ~isempty(savestr) ;

      if iscreenwrite == 1 ;
        disp('    BatchAdjustSM saving results') ;
      end ;

      orient landscape ;
      print('-dpsc2',[savestr 'DWDrepeatProj']) ;

      if iscreenwrite == 1 ;
        disp('    BatchAdjustSM finished save') ;
        disp('  ') ;
      end ;

    end ;


  end ;



  if viplot(4) ;    %  then make "after" 2-d scatterplot matrix

    if fignum == 1 ;    %  then this figure is first
      figure(1) ;
      clf ;
    else ;    %  then have an earlier figure, so shift this one
      vpos = get(gcf,'Position') ;
      vpos = [(vpos(1) + npixshift) ...
              (vpos(2) - npixshift) ...
               vpos(3) vpos(4)] ;
      figure(fignum) ;
      clf ;
      set(gcf,'Position',vpos) ;
    end ;
    set(gcf,'Name','PCA 2-d Projections, After') ;
    fignum = fignum + 1 ;

    mlegendcolorBatch = [[1 0 0] ; ...
                        [0 0 1]] ;

    mcolorBatch = [(batchlabels == 1)' zeros(n,2)]+...
                  [ zeros(n,2) (batchlabels == -1)'];


    paramstruct = struct('viout',[3], ...
                         'vipcplot',1:npc, ...
                         'icolor',mcolorBatch, ...
                         'iscreenwrite',iscreenwrite) ;

    if ~isempty(savestr) ;
      paramstruct = setfield(paramstruct, ...
                            'savestr',[savestr 'PCA2dAfter']) ;
    end ;

    if ~isempty(legcellstr) ;
      paramstruct = setfield(paramstruct, ...
                            'legendcellstr',legcellstr) ;
      paramstruct = setfield(paramstruct, ...
                            'mlegendcolor',mlegendcolorBatch) ;
    end ;

    curvdatSM(adjdata,paramstruct) ;


  end ;



end ;




