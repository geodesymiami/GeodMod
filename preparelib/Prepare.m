function [plotdataopt]=Prepare(opt)
%PreparePlot           - reads parameters from *min file and generates options for PlotData
%
%  usage:  [plotdataopt]=PreparePlot(opt)
%
%          'demfile'           Digtal Elevation Model                       [default 'off']
%          'ShadedRelieffile'  jpeg generated for now with demstuff_N       [default 'off']
%          'gpsfile'           gpsdata                                      [default 'off']
%          'quakefile'         earthquake location file                     [default 'off']
%          'quakefileformat'   earthquake file format ('anss',..)           [default 'anss']
%          'linefile'          faults                                       [default 'off']
%
%          NOTE: .rsc file are required for the dem and shaded relief *hgt.rsc and *jpeg.rsc
%          TODO: cleaner programming by reading as default the *dem file (2 byte integer) in *template
%          TODO: (Dem2Basemap could check from the filesize whether the input is a *dem or *hgt file)
%          TODO: notsure that the quakefile format is indeed anss. Need to find out name of format used.
%                        
%  Part of the MakeDataModel  suite
%  V1.0  Falk Amelung, September 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultopt=struct(                                                         ...
        'demfile'            ,        'off'                   ,            ...
        'ShadedRelieffile'   ,        'off'                   ,            ...
        'DemAreaEqualsSubset',        'off'                   ,            ...
        'subset'             ,        'off'                   ,            ...
        'CLim'               ,        'off'                   ,            ...
        'Cmap'               ,        'off'                   ,            ...
        'Fringe'             ,        'off'                   ,            ...
        'colorbaropt'        ,        'off'                   ,            ...
        'symbopt'            ,        'off'                   ,            ...
        'quakefile'          ,        'off'                   ,            ...
        'focalfile'          ,        'off'                   ,            ...
        'faultfile'          ,        'off'                   ,            ...
        'quakefileformat'    ,        'anss_readable'                  ,            ...
        'focalfileformat'    ,        'Focal'                          ,            ...  
        'ShadeFac'           ,        'off'                   ,            ...
        'GPStxtfile'         ,        'off'                   ,            ...
        'GPSfile'            ,        'off'                   ,            ...
        'GPSformat'          ,        'BenBrooks'             ,            ...
        'linefile'           ,        'off'                   ,            ...
        'faults_shapefile'   ,        'off'                   ,            ...
        'faults_shapefile2'  ,        'off'                   ,            ...
        'roads_shapefile'    ,        'off'                   ,            ...
        'GPSsta'             ,        'off'                   ,            ...
        'Sites'              ,        'off'                   ,            ...
        'Profile'            ,        'off'                   ,            ...
        'Quakes'             ,        'off'                   ,            ...
        'Focals'             ,        'off'                   ,            ...
        'Faults'             ,        'off'                   ,            ...
        'Roads'              ,        'off'                   ,            ...
        'ShadeOnly'          ,        'off'                   ,            ...
        'x_unit'             ,         'degrees'              ,            ...
        'marker1D'           ,         '{}'                     ,            ...
        'marker2D'           ,         '{}'                     ,            ...
        'FlipScale'          ,        'off'      )            ;
global dir_out 
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],defaultopt); end ;
[opt]=process_defaultoptions(opt,defaultopt);  %display(opt)
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end

%if  CheckInOut('','')  return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

makebasemapopt.demfile             = demfile;
makebasemapopt.ShadedRelieffile    = ShadedRelieffile;
makebasemapopt.DemAreaEqualsSubset = DemAreaEqualsSubset;
makebasemapopt.subset              = subset;
makebasemapopt.out_name            = fullfile(dir_out,'basemap.mat');

[new_flag]=Dem2Basemap(makebasemapopt);

load(makebasemapopt.out_name);
opt.basemap = basemap ;

% 
% read fault data
%
if linefile                                   % this is only kept for compability with existing Hawaii linefile
   load(linefile) ;
   Faults      = struct;
   Faults.lola = Lllh;
   Faults.xy   = lola2xy(Faults.lola',basemap,1)';
end

if faults_shapefile
       %Faults = MakeFaults(faultfile,basemap);  
       tmpfaults   = shaperead(faults_shapefile);
       Faults      = struct;
       Faults.lola = [[tmpfaults.X]' [tmpfaults.Y]'];
       Faults.xy   = lola2xy(Faults.lola',basemap,1)';
          if faults_shapefile2
             tmpfaults   = shaperead(faults_shapefile2);
             Faults.lola = [ Faults.lola ; [[tmpfaults.X]' [tmpfaults.Y]'] ];
             Faults.xy   = [  lola2xy(Faults.lola',basemap,1)'];
          end
end

%
% read roads shapefile
%
if roads_shapefile
       tmproads = shaperead(roads_shapefile);
       Roads       = struct;
       Roads.lola  = [[tmproads.X]' [tmproads.Y]'];
       Roads.xy    = lola2xy(Roads.lola',basemap,1)';        
end

% 
% read gps data
%
if GPStxtfile(1) || GPSfile(1)
  [fname,GPS]  = MakeGPS(opt);
  GPSsta = [GPS.lola;  GPS.enu]';
end

% 
% read quake data
%
if quakefile
    Quakes = MakeQuakes(quakefile,quakefileformat,basemap);   
end

%%%%%%% Introducing Focal Mechanism %%%%%%%%%%%%
if focalfile
	Focals  = MakeFocals(focalfile,focalfileformat,basemap);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isstruct(Profile)
   if ~isfield(Profile,'lola') && isfield(Profile,'xy')   Profile.lola=lola2xy(Profile.xy',basemap,-1)';  end
   if ~isfield(Profile,'xy') && isfield(Profile,'lola')   Profile.xy=lola2xy(Profile.lola',basemap, 1)';  end
   %if ~isfield(Profile,'xy') && isfield(Profile,'lalostrike')  
   %                                                      tmplola   =
   %                                                      Profile.xy=lola2xy(Profile.lola',basemap, 1)';  end
   %Profile.theta   = rad2deg(atan2(Profile.xy(2,2)-Profile.xy(1,2),Profile.xy(2,1)-Profile.xy(1,1)));     % ccw from horizontal
   vec              = [ Profile.xy(2,1)-Profile.xy(1,1),Profile.xy(2,2)-Profile.xy(1,2) ] ;     % ccw from horizontal
   Profile.unitvec  =  vec / norm(vec);
end

[marker1D,marker2D] = initialize_markers(marker1D,marker2D,basemap);

%TODO: the following probably works simeply with 'if quakefile'  and 'if linefile'
                          plotdataopt.Dem             = basemap.data   ;
                          plotdataopt.Profile         = Profile        ;
                          plotdataopt.Dem_unit        = basemap.Unit   ;
                          plotdataopt.x_unit          = x_unit         ;
                          plotdataopt.basemap         = basemap        ;         %TODO: May 2007: all code should use plotdataopt.basemap.x_first, etc 
if ShadeOnly              plotdataopt.ShadeOnly       = ShadeOnly      ;   end   % need Faults(1) in the case it is NaN

if GPSsta                 plotdataopt.Vectors         = GPSsta         ;   end
if Sites                  plotdataopt.Sites           = Sites          ;   end
if CLim                   plotdataopt.CLim            = CLim           ;   end
if Cmap                   plotdataopt.Cmap            = Cmap           ;   end   % Oct 1 2006: Did not work
if Fringe                 plotdataopt.Fringe          = Fringe         ;   end   % Oct 1 2006: Did not work
if ShadeFac               plotdataopt.ShadeFac        = ShadeFac       ;   end
if FlipScale              plotdataopt.FlipScale       = FlipScale      ;   end
if isstruct(Profile)      plotadataopt.Profile        = Profile        ;   end
if isstruct(Faults)       plotdataopt.Faults          = Faults         ;   end   % need Faults(1) in the case it is NaN
if isstruct(Roads)        plotdataopt.Roads           = Roads          ;   end
if isstruct(Quakes)       plotdataopt.Quakes          = Quakes         ;   end        
if isstruct(Focals)       plotdataopt.Focals          = Focals         ;   end
if isstruct(symbopt)      plotdataopt.Symbols         = symbopt        ;   end        
if isstruct(colorbaropt)  plotdataopt.colorbaropt     = colorbaropt    ;   end        
if iscell(marker1D)       plotdataopt.marker1D        = marker1D       ;   end   
if iscell(marker2D)       plotdataopt.marker2D        = marker2D       ;   end        
%
%  Generate baseplot
%
baseplotopt=plotdataopt;

out_name                     = fullfile(dir_out,'baseplot');
baseplotopt.colorbaropt.Title= 'Topography';
baseplotopt.Fringe           = false;

if CLim baseplotopt          = rmfield(baseplotopt,'CLim'); end
if Cmap baseplotopt          = rmfield(baseplotopt,'Cmap'); end

if new_flag
   logplot('PlotData',out_name,basemap,baseplotopt);
end
