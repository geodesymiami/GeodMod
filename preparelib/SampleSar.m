function []=SampleSar(opt)
%SampleSar                - generates quadtree'ed or grid dataset from rate
%
%  usage:  SampleSar(opt)
%
%          'dir_in'            Input directory (contains e.g. RsatA3)         [default 'off']
%          'igram_name'        filename for igram file                        [default is obtained with extract_name_fromSOdir]
%          'Medfiltsize'       median filtering prior to Quadtree/Gridding    [default 5]
%          'QtFittype'         Fittype (0 avg, 1  bilinear)(not testest)      [default 0]
%          'QtStartlevel'      Largest  quadrant has length 512/2^Startlevel  [default 2]
%          'QtEndlevel'        Smallest quadrant has length 512/2^Endlevel    [default 7]
%          'GridRowsCols'      [rows cols] for gridding                       [default 25]
%          'CoordSystem'       disloc: flipupdown and km for cx,coord,posting [default 'disloc']
%          'Plot'              plot data ? (screen,jpeg,ps,off)               [default 'off']
%
%  Sep 2006: I am not sure that QtEndlevel works well for anything else than 13 (default in original program). In some
%            For small QtEndlevel quadtree_part returns fields with sqval=10 which is the initial value hardwired in
%            quadtree_part. I do not know what this means. 
%  Attention: needs the capability to add a mask and remove all points located within this mask
%  For a grid size of 512, QtStartlevel=QtEndlevel=5 would produce a 32*32 grid each with 16 pixel
%
%  Part of the TimeSeries suite
%  V1.0  Falk Amelung, September 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global plot_save;
defaultopt=struct(                                                ...
        'DoIt'                       ,  'on'         ,            ...
        'in_name'                    ,  'off'        ,            ...
        'out_name'                   ,  'off'        ,            ...
        'Method'                     ,  'Quadtree'   ,            ...
        'Medfiltsize'                ,   1           ,            ...
        'QtTolerance'                ,   0.002       ,            ...
        'QtStartlevel'               ,   4           ,            ...
        'QtEndlevel'                 ,   7           ,            ...
        'QtFittype'                  ,   0           ,            ...
        'GridRowsCols'               ,   25          ,            ...
        'manualSampledFile'          ,   0           ,            ...
        'manualSampledFormat'        ,  'geodmod'    ,            ...
        'manualSampledPlotSaveFlag'  ,  'off'        ,            ...
        'plotdataopt'                ,  'off'        ,            ...
        'Plot'                       ,  'off'        )            ;
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],defaultopt); end ;
[opt]=process_defaultoptions(opt,defaultopt);           %display(opt)
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end
if  ~DoIt                         return; end
if  CheckInOut(in_name,out_name)  return; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(in_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[junk,junk,DataSet] = extract_ProjectName(in_name) ;
[radarlook]         = extract_hardwired_satparameters(in_name,'LOSvector') ;
if isfield(motion,'incangle')
   incangle            = motion.incangle;
   azimuth             = motion.azimuth;
else
   incangle          = extract_hardwired_satparameters(in_name,'incangle');
   azimuth           = extract_hardwired_satparameters(in_name,'azimuth');
   incangle          = ones(size(motion(1).data))*incangle;
   azimuth           = azimuth;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% median filtering of data
%
% FA March 2008. I believe the following can go because nanmedfilt2 is done in Igram2Motion. we do need tmpdata=motion.data;
tmpdata                = motion.data;
tmpamp                 = ones(size(tmpdata));
tmpamp(isnan(tmpdata)) = nan;

tmpdata(isnan(tmpamp)) = 10000;
tmpdata                = medfilt2(tmpdata,[Medfiltsize Medfiltsize]);
tmpdata(tmpdata==10000)= NaN;
motion.data            = tmpdata;     % Falk Jun 2007.  Not sure whether needed. Added because ProjectProfile introduced.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch Method
case{'Quadtree'}
   %
   % quadtree partioning
   %
   [newindmat,sqval,cx,cy,cntp,matsize] = quadtree_part(tmpdata,QtTolerance,QtFittype,QtStartlevel,QtEndlevel);

   nProb=find(sqval==10);  
   if nProb > 0
      logmessage(sprintf('Quadtree returns %d field with intial value 10. This may be a problem',length(nProb)))
      sqval(nProb)=NaN;               % SETTING SOME FIELDS TO NAN. NOT SURE THIS IS CORRECT (Sep 2006)
   end
   
   nN = find(~isnan(sqval));

   indmatnN = newindmat(nN,:); cntpnN = cntp(nN,:); cxnN = cx(:,nN); cynN = cy(:,nN); sqvalnN=sqval(nN);
   figure,subplot(2,2,1) ; imagesc(tmpdata); axis ij ; axis image; subplot(2,2,2) ; patch(cxnN,cynN,sqvalnN); axis ij ; axis image;  

   % here the unit is assigned to coord and cx,cy
   coord=[cntpnN(:,2)  cntpnN(:,1)]';
   cx=cxnN; cy=cynN;
   Ndata=length(sqvalnN);
   datavec=sqvalnN;

   % extract incangle at square centers
   tmpincangle = blkdiag(incangle,zeros(matsize-size(tmpdata,1),matsize-size(tmpdata,2)));
   incanglevec = zeros(Ndata,1);
   for i=1:Ndata;
        incanglevec(i) = nansum(nansum(tmpincangle(cy(2,i)+1:cy(3,i),cx(1,i)+1:cx(2,i))))/sum(sum(isfinite(tmpincangle(cy(2,i)+1:cy(3,i),cx(1,i)+1:cx(2,i)))));
       % there is a slight difference between the value obtained this way and from quadtree_part. Need to be investigated. For now we simply
       % use this as it is only for the incangle
   end

   logmessage(sprintf('Tolerance, Ndata  %6.4f %5d \n',QtTolerance,Ndata))        

   %%% clean data with the 2 commented lines %%%%%%%%%%
   %  patch(cx,cy,datavec); axis xy ; axis image;  
   % qq=deletepatch(cx,cy,datavec); datavec=qq;
   nN = find(isnan(datavec)==0); coord = coord(:,nN); cx = cx(:,nN); cy = cy(:,nN); datavec=datavec(nN); Ndata=length(datavec);
   %patch(cx,cy,datavec); axis xy ; axis image ;  
   % FA 1/2011: We may want to use Sjonni's interactive removepatches.m to clean the sampled data. We should load masks for cleaing as well.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case{'Grid'}
  if length(GridRowsCols) == 1
     GridRows=GridRowsCols;
     GridCols=GridRowsCols;
  elseif length(GridRowsCols) == 2
     GridRows=GridRowsCols(1);
     GridCols=GridRowsCols(2);
  else
    error('GridRowsCols')
  end

  gridsize_row =size(tmpdata,1)/GridRows ;   
  gridsize_col =size(tmpdata,2)/GridCols ;   

  yvec=[0:GridRows-1];
  xvec=[0:GridCols-1];
  [yind,xind]=deal(zeros(1,GridRows*GridCols));
  [yind(:),xind(:)]=meshgrid(yvec,xvec);
  yind=yind*gridsize_row + gridsize_row/2 ;
  xind=xind*gridsize_col + gridsize_col/2 ;

  cy=[yind-gridsize_row/2; yind-gridsize_row/2; yind+gridsize_row/2; yind+gridsize_row/2] ;
  cx=[xind-gridsize_col/2; xind+gridsize_col/2; xind+gridsize_col/2; xind-gridsize_col/2] ;
  
  datavec     = tmpdata(sub2ind(size(tmpdata),round(yind),round(xind)));
  incanglevec = incangle(sub2ind(size(tmpdata),round(yind),round(xind))); 
 
  % remove not-unwrapped points
  ind_nan              = find(isnan(datavec)) ;
  datavec(ind_nan)     = [];
  incanglevec(ind_nan) = [];
  
  yind(ind_nan)        = [];
  xind(ind_nan)        = [];
  cx(:,ind_nan)        = [];
  cy(:,ind_nan)        = [];

  Ndata=length(datavec);
  coord=[xind ; yind];
  cx   = cx ;
  cy   = cy ;
    
case{'ProjectProfile'}
  profiledata           = MakeProfile(motion,profileopt);
  %[projvec,profiledata] = MakeProfile(motion,profileopt);
  %radarlook       = radarlook * projvec ;

  coord            = zeros(2,size(profiledata.profileall,1));
  coord(1,:)       = profiledata.profileall(:,1) ;
  coord(2,:)       = repmat(profileopt.projectpara.dir,1,length(coord)) ;
  datavec          = profiledata.profileall(:,2)';
 
  cx               = [];                                      % ToDo: If MakeProfile returns cx,cy we could plot using patch below 
  cy               = [];
  Ndata            = length(datavec); 
  
  logplot('PlotProfile',out_name,coord(1,:),datavec)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Assign output fields
%

amp = ones(size(motion.data)); amp(isnan(motion.data)) = nan;

dataset.data        = motion.data ;
dataset.amp         = amp  ;                  % not sure whether really needed. Try to remove at later stage.
dataset.x_first     = motion.x_first  ; 
dataset.y_first     = motion.y_first  ; 
dataset.x_step      = motion.x_step  ; 
dataset.y_step      = motion.y_step  ; 
dataset.x_unit      = motion.x_unit  ; 
dataset.Unit        = motion.Unit    ; 
dataset.TotalTime   = motion.TotalTime;
dataset.StartDate   = motion.StartDate;
dataset.radarlook   = repmat(radarlook',1,Ndata);
dataset.radarlook   = GenerateLOSVec(incanglevec,azimuth-90.0);
dataset.Ndata       = Ndata;
dataset.datavec     = datavec ;
dataset.coord       = coord ;  if strcmp(Method,'ProjectProfile') dataset.CoordSystem = Method ;  end
dataset.cx          = cx   ;
dataset.cy          = cy   ;
dataset.DataSet     = DataSet ;
dataset.tolerance   = QtTolerance ;

if manualSampledFile
  if  strmatch(manualSampledFormat,'Sjonni')
      S             = load(manualSampledFile);
      g             = fieldnames(S);
      
      dataset       = qtsjonni2qtgeodmod(S.(g{1}), dataset, plotdataopt.basemap);
      
      %dataset.coord = S.(g{1}).cnt;
      %dataset.datavec= S.(g{1}).sv';
      %dataset.cx    = S.(g{1}).cx;
      %dataset.cy    = S.(g{1}).cy;
      %dataset.radarlook= S.(g{1}).los;
      %dataset.cov      = S.(g{1}).cov;
      %dataset.Ndata    = length(dataset.datavec);

  else
      load(manualSampledFile);       % Read fixed Qt file in geodmod format
  end
end
% ToDo: call patch from a higher-level program such as PlotPatch that generates the figure handle, axis to simplyply logplot
%patch(cx,cy,datavec); axis xy ; axis image ;  
if sum(strcmp(Method,{'Quadtree','Grid'}))
    tmpplotdataopt           = plotdataopt;
    tmpplotdataopt.PlotType  = 'SampledData';
    tmpplotdataopt.ShadeOnly = true;
    
    tmpplot_save             = plot_save;
    if ~manualSampledPlotSaveFlag plot_save=false; end       % FA 1/2011: Switching off plotting because it is very slow
    logplot('PlotData',out_name,dataset,tmpplotdataopt)      % FA 9/2008  removed because it crashed matlab. Don't know why. Seems to be related to the plot status matlab is in.
    plot_save                = tmpplot_save;
end
save(out_name,'dataset') ;
