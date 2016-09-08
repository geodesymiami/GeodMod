function [dataset]=MakeDataset(opt)
%MakeDataset             - make dataset for modelling from  quadtree'ed or gridded interferograms and GPS data
%
%  usage:  [dataset]=MakeDataset(opt)
%
%          'insarlist'         One or more quadtree'ed datasets (*Qt*mat)   
%          'GPSfile'           GPS data set                                 [default 'off']
%          'GPShorz'           Add horizontal GPS velocities to dataset     [default 'off']
%          'GPSvert'           Add vertical GPS velocities to dataset       [default 'off']
%                              (only used if no SAR data given)             
%          'SARsigphi'         uncertainty for SAR (1,'SARequal',equal)     [default 1]
%          'GPSsigphi'         uncertainty for GPS                          [default 1]
%          'Plot'              plot data ? (screen,jpeg,ps,off)             [default 'off']
%                        
%  Part of the TimeSeries suite
%  V1.0  Falk Amelung, September 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  set default options. Process options fom input or from file *.min
defaultopt=struct(                                                         ...
        'DoIt'               ,        'on'                    ,            ...
        'insarlist'          ,        'off'                   ,            ...
        'GPSfile'            ,        'off'                   ,            ...
        'GPShorz'            ,        'off'                   ,            ...
        'GPSvert'            ,        'off'                   ,            ...
        'Unit'               ,        'm/yr'                  ,            ...
        'sigphi'             ,        1                       ,            ...
        'SARsigphi'          ,        1                       ,            ...
        'GPSsigphi'          ,        1                       ,            ...
         'CoordSystem'        ,        'Disloc'               ,            ...
        'Plot'               ,        'off'      )            ;

if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],defaultopt); end ;
[opt] = process_defaultoptions(opt,defaultopt);  %display(opt)
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end

if ~DoIt   dataset=[];      return; end
if  CheckInOut(insarlist,'')  return; end

if ~iscell(insarlist)  insarlist=cellstr(insarlist);  end    % convert insarlist into cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load SAR data and fill into dataset structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if length(insarlist)==0 tmpdataset=[]; end
 for i=1:length(insarlist) 
     insarfile=insarlist{i};
     load(insarfile) ; tmpdataset(i)=dataset;
 end
 dataset=tmpdataset ; clear tmpdataset ;

  x_posting = basemap.x_posting ;
  y_posting = basemap.y_posting ;
  
  if isfield(dataset(i),'CoordSystem') CoordSystem = dataset(i).CoordSystem ;  end
  
  for i=1:length(dataset)
      dataset(i).x_posting = x_posting  ;
      dataset(i).y_posting = y_posting  ;
      dataset(i).CoordSystem = CoordSystem ;
      dataset(i).datavec    =double(dataset(i).datavec);
  end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% flipud data for Disloc 
%
switch CoordSystem
    case{'Disloc'}
           for i=1:length(dataset)              % FA Feb 2008: if no SAR data length is zero and this step is skipped
               dataset(i).coord(2,:)  = size(dataset(i).data,1) - dataset(i).coord(2,:) ;
               dataset(i).cy          = size(dataset(i).data,1) - dataset(i).cy         ;
           end
       %
       % convert coord,cx,cy into km. flipud data (June 2007: This is probably only for PlotModel,not sure)
       %
           for i=1:length(dataset)
               dataset(i).coord     = [dataset(i).coord(1,:)*x_posting; dataset(i).coord(2,:)*y_posting ] ;  %mogi and disloc require EN coordinates
               dataset(i).cx        = dataset(i).cx * x_posting ;
               dataset(i).cy        = dataset(i).cy * y_posting ;
            end
    case{'ProjectProfile'}
       % nothing to be done becuse coord is already in km
    otherwise
       errordlg('CoordSystem: %s not allowed, -- exiting',CoordSystem);
end

% TODO: Cleaner programming would be to do the Discloc CoordSystem conversion for
%       everything (SAR,GPS,DEM) in one block at the end of this function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load GPS file and add to dataset structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if GPSfile(1)
   load(GPSfile)
   if GPShorz
      i=length(dataset)+1;
      [dataset(i).DataSet      ,dataset(i+1).DataSet      ] = deal( 'GPSeast','GPSnorth'                  );
      [dataset(i).datavec      ,dataset(i+1).datavec      ] = deal( [GPS.e_rate],  [GPS.n_rate]           );  %TODO: Change e_rate->e
      [dataset(i).normalization,dataset(i+1).normalization] = deal( [GPS.e_error], [GPS.n_error]          );  %TODO: Change e_rate->e
      [dataset(i).Ndata        ,dataset(i+1).Ndata        ] = deal( length(dataset(i).datavec)            );
      [dataset(i).coord        ,dataset(i+1).coord        ] = deal( [ [GPS.x]; [GPS.y] ]                  );
      [dataset(i).dcov         ,dataset(i+1).dcov         ] = deal( GPS(1).cov                            );
      [dataset(i).Unit         ,dataset(i+1).Unit         ] = deal( GPS(1).Unit                           );
      [dataset(i).CoordSystem  ,dataset(i+1).CoordSystem  ] = deal( CoordSystem                           );
      
      dataset(i).radarlook                                  =  repmat([1 0 0]',1,length(dataset(i).datavec));
      dataset(i+1).radarlook                                =  repmat([0 1 0]',1,length(dataset(i).datavec));
      
  end
   if GPSvert
      i=length(dataset)+1;
      [dataset(i).DataSet                                 ] = deal( 'GPSvert'                             );
      [dataset(i).datavec                                 ] = deal( [GPS.u_rate]                          ); %TODO: Change e_rate->e
      [dataset(i).normalization                           ] = deal( [GPS.u_error]                         ); %TODO: Change e_rate->e
      [dataset(i).Ndata                                   ] = deal( length(dataset(i).datavec)            );
      [dataset(i).coord                                   ] = deal(    [ [GPS.x]; [GPS.y] ]               );

      [dataset(i).dcov                                    ] = deal( GPS(1).cov                            );
      [dataset(i).Unit                                    ] = deal( GPS(1).Unit                           );
      [dataset(i).CoordSystem                             ] = deal( CoordSystem                           );
      
       dataset(i).radarlook                                 = repmat([0 0 1]',1,length(dataset(i).datavec));
   
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add elevation to dataset(for inversion using topo approximation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch CoordSystem
       case{'Disloc' 'ProjectProfile'}
           dem   = flipud(basemap.data);
       otherwise
           error('CoordSystem: %s not allowed, -- exiting',CoordSystem);
end

if strcmp(basemap.Unit,'Meters')
      dem   = dem/1000;                  % convert because source depth is in km
   else
       error('user error: Unit %s not recognized -- exiting', basemap.Unit) ;
end

if exist('x_posting','var')
    [dataset]=extract_hgt_from_dem(dem,dataset,x_posting,y_posting); 
else
    [dataset]=extract_hgt_from_dem(dem,dataset); 
end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine whether SAR  or GPS,
% assign dataset.SAR, dataset.GPS and dataset.ind  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
last=0;
for i=1:length(dataset)
    dataset(i).SAR= isempty(strmatch('GPS',dataset(i).DataSet));
    dataset(i).GPS=~isempty(strmatch('GPS',dataset(i).DataSet));
    dataset(i).ind=last+1:last+dataset(i).Ndata;
    dataset(i).exist=true;
    last=last+dataset(i).Ndata;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add sigphi to dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(dataset)
    if isstruct(sigphi)
       if dataset(i).SAR
           dataset(i).sigphi = sigphi.SAR{i};
       else
           DataSetStr=dataset(i).DataSet;if strcmp(DataSetStr,'GPSeast') | strcmp(DataSetStr,'GPSnorth') DataSetStr='GPShorz';end
           dataset(i).sigphi=sigphi.(DataSetStr);        % FA 1/11: for GPShorz or GPSvert
       end
    else
       dataset(i).sigphi=sigphi;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add Ramp to dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('ramp') 
   if ~(length(ramp)==sum([dataset.SAR])) error('Need ramp parameter for all SAR data sets -- existing'); end
   m_phaseramp = [ramp{:}];
        tmpobjfuncopt.PhaseRamp = true;
        tmpobjfuncopt.FactorLin = false;   
        [dataset]=ModifyDatasetLin(dataset,tmpobjfuncopt,m_phaseramp');
        %G_phaseramp_coord = [ ones(dataset(i).Ndata,1)            dataset(i).coord(1:2,:)'];
end
