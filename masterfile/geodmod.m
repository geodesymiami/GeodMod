function []=geodmod(fname,varargin)
% GEODMOD  -  reads and save roi_pac-processed interferograms,gps data and models data
%
%  usage:  geodmod fname.min [commands] [opt_commands]
%          geodmod fname
%
%  e.g.    geodmod Darwin.min
%          geodmod Darwin
%          geodmod Darwin gclean model 
%          geodmod Darwin gclean model close all 
%          geodmod Darwin plotsurface3dopt.z_offset=1
%          geodmod Darwin dir_out_parent='''.'''                % to run distribution examples locally     
%          geodmod Darwin gclean model close all modelopt.N_mogi=2 modelopt.N_disloc=1
%
%          need to run addpath_geodmod  to set the directories
%          
%  to learn about options:
%     help Prepare
%     help MakeInSAR
%     help MakeGPS
%     help InverseNonLin
%     help InverseLin
%
%     help gclean
%     
%     e.g. gclean model
%          gclean all
%
%     gclean remembers the name of the last output directory (dir_out in .last_dir_out)
%
%    Quick testing by running:
%
%    geodmod Darwin_classexample.min
%    geodmod Forward_mogi.min
%    geodmod MLGPSonly.min
%    geodmod basemap_hawaii.min
%    geodmod ML2002-2005.min
%    test_geodmod
%    test_geodmod 'all'
%  
% June 2007 Falk Amelung

  clear  dir_out plot_visible plot_save 
  global dir_out plot_visible plot_save 
  
 if nargin==0 help geodmod; return; end 
 
%
% add .min to fname if necessary
%
 
 [path,name,ext] = fileparts(fname);
 if isempty(ext) fname=fullfile(path,[name '.min']); end
 dir_out         = name;                                    % set out_dir to filename [default]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% Initialize some variables,  read *min file %
%                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 set(0,'DefaultPatchLineStyle','none')
 plot_visible               = true; 
 plot_save                  = 'pdf'; 
 makesaropt.DoIt            = 'off';
 makegpsopt.DoIt            = 'off';
 inverseopt.DoIt            = 'off';
 plotthemodelopt            = [];
 plotsurface3dopt           = [];
 makesaropt.data2igramopt   = [];
 modelopt                   = [];

 [opt]=read_options_fromfile(fname,[]); f=fieldnames(opt); for i=1:length(f); eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end ; 

 for i=1:length(varargin)                                    % manual options given on command line
      if strfind(varargin{i},'=')
           C       = textscan(strtrim(varargin{i}),'%s%s','delimiter','=');
           keyword = char(C{1});
           value   = char(C{2});
           if any(isalpha(value)) && ~any(strfind(value,''''))
               value_mod = sprintf('''%s''',value);
           else
               value_mod = value;
           end
           eval(sprintf('%s=%s;',keyword,value_mod));
      end       
 end
 
 if ~isfield(makesaropt,'plotdataopt');      makesaropt.plotdataopt     =[]; end
 if ~isfield(inverseopt,'plotdataopt');      inverseopt.plotdataopt     =[]; end
 if ~isfield(inverseopt,'plotmodelopt');     inverseopt.plotmodelopt    =[]; end                       % FA 7/08: needed as long PlotModel is used
 if ~isfield(plotthemodelopt,'plotdataopt'); plotthemodelopt.plotdataopt=[]; end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                             %             
% Determine out_dir, do .last_dir_out, gclean %
%                                             %             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if exist('dir_out_parent','var') && ~isempty(dir_out_parent) 
    dir_out=[dir_out_parent filesep dir_out] ;                     % use different dir_out_parent if given
 end   

 logmessage(sprintf('output directory: %s ',dir_out))
 f1=fopen('.last_dir_out','w'); fprintf(f1,dir_out); fclose(f1);

 i_gclean = strmatch('gclean',varargin);
 i_close  = strmatch('close',varargin);
 
 tmpvarargin = varargin;
 if i_gclean 
     logmessage('cleaning...' );
     gclean(tmpvarargin{i_gclean+1});
     tmpvarargin(i_gclean:i_gclean+1)='';
 end

 if i_close
     close(tmpvarargin{i_close + 1});
     tmpvarargin(i_close:i_close+1)='';
 end
 
if ~exist(dir_out,'dir')  mkdir(dir_out) ; end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% Check for Data Dirs, Download if necessary %
%                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %CheckDownloadData(prepareopt,makesaropt) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% START DOING THINGS %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% read DEM and prepare plotdataopt           %
%                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
 % select the same subset as from the interferogram  
 if exist('makesaropt','var') && isfield(makesaropt,'subset');  prepareopt.subset=makesaropt.subset; end

 [plotdataopt]              = Prepare(prepareopt);
 makesaropt.plotdataopt     = process_defaultoptions(makesaropt.plotdataopt,plotdataopt);
 inverseopt.plotdataopt     = process_defaultoptions(inverseopt.plotdataopt,plotdataopt);
 plotthemodelopt.plotdataopt= process_defaultoptions(plotthemodelopt.plotdataopt,plotdataopt);
 makegpsopt.basemap         = plotdataopt.basemap;
 makedatasetopt.basemap     = plotdataopt.basemap;
 
 if  isfield(inverseopt,'distribopt'); inverseopt.plotdataopt.distribopt = inverseopt.distribopt; end                          % FA 9/08: makes PlotThresh is available in plotsurface3d
 if  isfield(inverseopt,'QuickStop') && strcmp(inverseopt.QuickStop,'on'); inverseopt.distribopt.DoIt=0; plotthemodelopt.DoIt=false;plotsurface3dopt.DoIt = false; end
 if  isfield(plotthemodelopt,'DoIt') && strcmp(plotthemodelopt.DoIt,'off'); plotsurface3dopt.DoIt = false; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% now make InSAR,GPS data. Prepare dataset.  %
%                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [makedatasetopt.insarlist]= MakeSAR(makesaropt);
 [makedatasetopt.GPSfile]  = MakeGPS(makegpsopt);
 [dataset]                 = MakeDataset(makedatasetopt);

 %if isfield (makedatasetopt,'GPSfile')   prepareopt.gpsfile=makedatasetopt.GPSfile; end     % This would plot GPS data on all PlotData plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% Set-up model, conduct inversions           %
%                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [modelopt] = InitializeModelopt(modelopt,plotdataopt.basemap);
 [modelopt] = InverseNonLin(dataset,inverseopt,modelopt);                                           
 [modelopt] = InverseLinDistrib(dataset,inverseopt,modelopt);                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% Plot model predictions (2D and 3D)         %
%                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 PlotTheModel(dataset,modelopt,plotthemodelopt);                                           
 PlotSurface3D(dataset,modelopt,inverseopt.plotdataopt,plotsurface3dopt);   

 disp(sprintf('Done: geodmod %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s',fname,varargin{:}))

% Other modelling software                      %%

 %ModelVisco1d(modelvisco1dopt);
