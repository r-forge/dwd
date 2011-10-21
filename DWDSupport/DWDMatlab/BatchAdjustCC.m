function adjdata = BatchAdjustCC(rawdata,batchlabels,paramstruct) 
% BATCHADJUSTCC, of adjustment of "centerpoint" of subpopulations,
%   Replaces Steve Marron's matlab function BatchAdjustSM, which adds 
%     graphics features
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
%    icolor           1  (default)  Red for Batch label = 1, Blue for Batch label = -1
%                     nx3 color matrix:  a color label for each Batch label
%
%    markerstr        Used for figures 1 & 4 only
%                     Can be either a single string with symbol to use for marker,
%                         e.g. 'o', '.', '+', 'x'
%                         (see "help plot" for a full list)
%                     (Default) '+' for Batch label = 1, 'o' for Batch label = -1
%                     Or a character array (n x 1), of these symbols,
%                         One for each data vector, created using:  strvcat
%
%    isubpopkde       0  construct kde using only the full data set
%                     1  (default) partition data into subpopulations, using the color
%                            indicators in icolor (defaults to 0, unless icolor
%                            is an nx3 color matrix), as markers of subsets.
%                            The corresponding mixture colors are then used in
%                            the subdensity plot, and overlaid with the full 
%                            density shown in black
%                     2  Show only the component densities (in corresponding 
%                            colors), without showing the full population
%                            density
%
%    idatovlay        0  Do not overlay data on kde plots (on diagonal)
%                     1  (default) overlay data using heights based on data ordering
%                              Note:  To see "c.d.f. style" increasing line, 
%                                     should also sort the data
%                     2  overlay data using random heights
%                     another integer > 0,  overlay data, using random heights,
%                                           with this numbers as the seed (so can 
%                                           better match data points across plots),
%                                           (should be an integer with <= 8 digits)
%
%    ndatovlay     number of data points overlayed (only has effect for idatovlay > 0)
%                       1  -  (default) overlay up to 1000 points 
%                                           (random choice, when more)
%                       2  -  overlay full data set
%                       n > 2   -  overlay n random points
%
%    datovlaymax      maximum (on [0,1] scale, with 0 at bottom, 1 at top of plot)
%                     of vertical range for overlaid data.  Default = 0.6
%
%    datovlaymin      minimum (on [0,1] scale, with 0 at bottom, 1 at top of plot)
%                     of vertical range for overlaid data.  Default = 0.5
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
%                          for Matlab 6.0, may be different for other
%                          Matlab versions)
%
%    mlegendcolor     2 x 3 color matrix, corresponding to cell legends above
%                     (defaults to red and blue when not specified)
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
%    npcadiradd       Number of Principal Component Directions to Add
%                     0  Don't add any, just directions in mdir
%                     integer > 0   Add this number of PC directions to mdir
%                         Note this assumes that no elements of mdir are PCA
%                         directions.  Otherwise will get an error.
%                     Default: 4
%                     integer < 0   Add this number of PC directions,
%                         for the projections on the subspace orthogonal
%                         to the subspace generated by mdir
%
%
% Output:
%   adjdata     - d x n matrix of adjusted data, 
%                          columns are cases, rows are genes 
%    

% Assumes path can find personal functions:
%    bwsjpiSM.m
%    kdeSM.m
%    vec2matSM.m
%    DWD1SM.m
%    curvdatSM.m
%    axisSM.m
%    sepelimdwd.m
%    lbinrSM.m
%    bwrfphSM.m
%    bwosSM.m
%    rootfSM.m
%    bwrotSM.m
%    bwsnrSM.m
%    iqrSM.m
%    cquantSM.m
%    pcaSM.m
%    sizerSM.m
%    rmeanSM.m
%    madSM.m
%    projplot1SM.m
%    qqLM.m
%    scatplotSM.m
%    bwrswSM.m
%    nprSM.m
%    sz2SM.m
%    sc2SM.m
%    CHkdeSM.m
%    CHsz1SM.m
%    CHlbinrSM.m
%    KMcdfSM.m
%    LBcdfSM.m


%    Copyright (c) J. S. Marron 2003, Chris Cabanski 2009-2010


%  First set paths to find needed subroutines
%
addpath subroutines -end ;
addpath subroutines\general -end ;
addpath subroutines\smoothing -end ;
addpath subroutines\BatchAdjust -end ;
addpath subroutines\SDPT3-4.0-beta\solver -end ;
addpath subroutines\SDPT3-4.0-beta\solver\Mexfun -end ;



%  Set all input parameters to defaults
%
viplot = [1; 1; 0; 1] ;
n = size(rawdata,2) ;
icolor = [(batchlabels == 1)' zeros(n,2)]+...
                  [ zeros(n,2) (batchlabels == -1)']; 
isubpopkde = 1 ;
idatovlay = 1 ;
ndatovlay = 1 ;
datovlaymax = 0.6 ;
datovlaymin = 0.5 ;
savestr = [] ;
titlestr = ['Batch Adjustment'] ;
titlefontsize = [] ;
legcellstr = {} ;
mlegendcolor = [[1 0 0]; [0 0 1]] ;
iscreenwrite = 0 ;
minproj = [] ;
maxproj = [] ;
npcadiradd = 4 ;
markerstr = [] ;
  for i = 1:length(batchlabels) ;
    if batchlabels(i) == 1 ;
      markerstr = strvcat(markerstr,'+') ;
    elseif batchlabels(i) == -1 ;
      markerstr = strvcat(markerstr,'o') ;
    else ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from BatchAdjustSM:          !!!') ;
      disp('!!!   batchlabels not in correct format    !!!') ;
      disp('!!!   i.e. batchlabels are not +1 and -1   !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    end ;
  end ;



%  Now update parameters as specified,
%  by parameter structure (if it is used)
%
if nargin > 2 ;   %  then paramstruct is an argument

  if isfield(paramstruct,'npcadiradd') ;    %  then change to input value
    npcadiradd = getfield(paramstruct,'npcadiradd') ; 
  end ;

  if isfield(paramstruct,'viplot') ;    %  then change to input value
    viplot = getfield(paramstruct,'viplot') ; 
  end ;

  if isfield(paramstruct,'icolor') ;    %  then change to input value
    icolor = getfield(paramstruct,'icolor') ; 
  end ;

  if isfield(paramstruct,'isubpopkde') ;    %  then change to input value
    isubpopkde = getfield(paramstruct,'isubpopkde') ; 
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

  if isfield(paramstruct,'markerstr') ;    %  then change to input value
    markerstr = getfield(paramstruct,'markerstr') ; 
  end ;

  if isfield(paramstruct,'idatovlay') ;    %  then change to input value
    idatovlay = getfield(paramstruct,'idatovlay') ; 
  end ;

  if isfield(paramstruct,'ndatovlay') ;    %  then change to input value
    ndatovlay = getfield(paramstruct,'ndatovlay') ; 
  end ;

    if isfield(paramstruct,'datovlaymax') ;    %  then change to input value
    datovlaymax = getfield(paramstruct,'datovlaymax') ; 
  end ;

  if isfield(paramstruct,'datovlaymin') ;    %  then change to input value
    datovlaymin = getfield(paramstruct,'datovlaymin') ; 
  end ;

if isfield(paramstruct,'mlegendcolor') ;    %  then change to input value
    mlegendcolor = getfield(paramstruct,'mlegendcolor') ; 
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
  np = sum(flagp) ;
  nn = sum(flagn) ;
  icolorp = icolor(flagp,:) ;
  icolorn = icolor(flagn,:) ;


  %  Find DWD direction
  %
  dirvec = DWD1SM(rawdata(:,flagp),rawdata(:,flagn)) ;


  %  Project data
  %
  vprojp = rawdata(:,flagp)' * dirvec ;
  vprojn = rawdata(:,flagn)' * dirvec ;

  meanprojp = mean(vprojp) ;
  meanprojn = mean(vprojn) ;


  %  Subtract respective class means
  %
  adjdata(:,flagp) = rawdata(:,flagp) - vec2matSM(meanprojp * dirvec,length(vprojp)) ;
  adjdata(:,flagn) = rawdata(:,flagn) - vec2matSM(meanprojn * dirvec,length(vprojn)) ;



  %  Make graphics (if needed)
  %
  fignum = 1 ;
  if viplot(1) ;    %  then make "before" 2-d scatterplot matrix

    figure(fignum) ;
    clf ;
    fignum = fignum + 1 ;
    set(gcf,'Name','PCA 2-d Projections, Before') ;

    mlegendcolorBatch = mlegendcolor ;


    paramstruct = struct('viout',[3], ...
                         'vipcplot',1:npcadiradd, ...
                         'isubpopkde',isubpopkde, ...
                         'idatovlay',idatovlay, ...
                         'ndatovlay',ndatovlay, ...
                         'datovlaymax',datovlaymax, ...
                         'datovlaymin',datovlaymin, ...
                         'markerstr',markerstr, ...
                         'icolor',icolor, ...
                         'mlegendcolor',mlegendcolorBatch, ...
                         'iscreenwrite',iscreenwrite) ;

    if ~isempty(savestr) ;
      paramstruct = setfield(paramstruct, ...
                            'savestr',[savestr 'PCA2dBefore']) ;
    end ;

    if ~isempty(legcellstr) ;
      paramstruct = setfield(paramstruct, ...
                            'legendcellstr',legcellstr) ;
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
    plot(xgrid1,vkdep1,'-','LineWidth',2,'Color',mlegendcolor(1,:)) ;
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
        plot(xgrid1,vkden1,'-','LineWidth',2,'Color',mlegendcolor(2,:)) ;

        plotbottom = vax(3) ;
        plottop = vax(4) ;
        yrand = plotbottom + (0.7 + 0.1 * rand(length(vprojp),1)) ...
                                              * (plottop - plotbottom) ;
               %  y coords for jittering
        for i = 1:np ;
           plot(vprojp(i),yrand(i),'+','Color',icolorp(i,:)) ;
        end;

        yrand = plotbottom + (0.5 + 0.1 * rand(length(vprojn),1)) ...
                                              * (plottop - plotbottom) ;
               %  y coords for jittering
        for i = 1:nn ;
           plot(vprojn(i),yrand(i),'o','Color',icolorn(i,:)) ;
        end;

        if ~isempty(legcellstr) ;
            tx = vax(1) + 0.1 * (vax(2) - vax(1)) ;
            ty = vax(3) + 0.9 * (vax(4) - vax(3)) ;
            voutstr = legcellstr{1} ;
          text(tx,ty,voutstr,'Color',mlegendcolor(1,:),'FontSize',18) ;
            ty = vax(3) + 0.4 * (vax(4) - vax(3)) ;
            vnotoutstr = legcellstr{2} ;
          text(tx,ty,vnotoutstr,'Color',mlegendcolor(2,:),'FontSize',18) ;
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
    plot(xgrid2,vkdep2,'-','LineWidth',2,'Color',mlegendcolor(1,:)) ;
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
        plot(xgrid2,vkden2,'-','LineWidth',2,'Color',mlegendcolor(2,:)) ;

        plotbottom = vax(3) ;
        plottop = vax(4) ;
        yrand = plotbottom + (0.7 + 0.1 * rand(length(vreprojp),1)) ...
                                              * (plottop - plotbottom) ;
               %  y coords for jittering
        for i = 1:np ;
           plot(vreprojp(i),yrand(i),'+','Color',icolorp(i,:)) ;
        end;

        yrand = plotbottom + (0.5 + 0.1 * rand(length(vreprojn),1)) ...
                                              * (plottop - plotbottom) ;
               %  y coords for jittering
        for i = 1:nn ; 
           plot(vreprojn,yrand,'o','Color',icolorn(i,:)) ;
        end;

        if ~isempty(legcellstr) ;
            tx = vax(1) + 0.1 * (vax(2) - vax(1)) ;
            ty = vax(3) + 0.9 * (vax(4) - vax(3)) ;
            voutstr = legcellstr{1} ;
          text(tx,ty,voutstr,'Color',mlegendcolor(1,:),'FontSize',18) ;
            ty = vax(3) + 0.4 * (vax(4) - vax(3)) ;
            vnotoutstr = legcellstr{2} ;
          text(tx,ty,vnotoutstr,'Color',mlegendcolor(2,:),'FontSize',18) ;
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

    mlegendcolorBatch = mlegendcolor ;


    paramstruct = struct('viout',[3], ...
                         'vipcplot',1:npcadiradd, ...
                         'isubpopkde',isubpopkde, ...
                         'markerstr',markerstr, ...
                         'idatovlay',idatovlay, ...
                         'ndatovlay',ndatovlay, ...
                         'datovlaymax',datovlaymax, ...
                         'datovlaymin',datovlaymin, ...
                         'icolor',icolor, ...
                         'mlegendcolor',mlegendcolorBatch, ...
                         'iscreenwrite',iscreenwrite) ;

    if ~isempty(savestr) ;
      paramstruct = setfield(paramstruct, ...
                            'savestr',[savestr 'PCA2dAfter']) ;
    end ;

    if ~isempty(legcellstr) ;
      paramstruct = setfield(paramstruct, ...
                            'legendcellstr',legcellstr) ;
    end ;

    curvdatSM(adjdata,paramstruct) ;


  end ;



end ;




diary
