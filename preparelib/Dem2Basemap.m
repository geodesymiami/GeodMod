function [new_flag]=Dem2Basemap(opt)
%Dem2Basemap        - reads dem and shaded relief *jpg  and makes basemap file
%
%  usage:  []=Dem2Basemap(opt)
%
%          'demfile'            Digtal Elevation Model                       [default 'off']
%          'ShadedRelieffile'   jpeg generated for now with demstuff_N       [default 'off']
%
%          NOTE: .rsc file are required for the dem and shaded relief *hgt.rsc and *jpeg.rsc
%          TODO: cleaner programming by reading as default the *dem file (2 byte integer) in *template
%          TODO: (Dem2Basemap could check from the filesize whether the input is a *dem or *hgt file)
%                        
%  Part of the MakeDataModel  suite
%  V1.0  Falk Amelung, September 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dir_out

defaultopt=struct(                                                        ...
        'demfile'            ,        'off'                  ,            ...
        'ShadedRelieffile'   ,        'off'                  ,            ...
        'subset'             ,        'off'                               ...
                 )          ;


if ~exist('opt','var');  [opt]=read_options_fromfile([mfilename '.min'],defaultopt); end 
[opt]=process_defaultoptions(opt,defaultopt);  %display(opt)
f=fieldnames(opt) ; for i=1:length(f); eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end

[pathstr, name, ext] = fileparts(ShadedRelieffile); 

inlist  = {demfile [demfile '.rsc'] ShadedRelieffile [ShadedRelieffile '.rsc'] };
outfile = fullfile(dir_out,'basemap.mat');

new_flag=true;
if  CheckInOut(inlist,outfile); new_flag=false;  return; end

%
% load DEM
%
   lopt.Unit = 'm';               % no sign changes
   lopt.DataType='Dem';
   lopt.origin='dem';
   if isstruct(subset); lopt.subset=subset ; elseif subset;  lopt.subset=subset;  end
   if DemAreaEqualsSubset 
      lopt=rmfield(lopt,'subset');
   end
 
   [basemap] = LoadData(demfile,lopt);

   [dem]=add_shade2Data(basemap,ShadedRelieffile,lopt);  %TODO: This function should be call load_shadedrelief_and_add2igram

   basemap.shade=dem.shade ;
   basemap.DataSet='Topography';
   basemap.Unit='Meters';      % TODO: We should have this in the rsc file for the dem and read the unit from there

%
%  Calculate x_posting,y_posting. Will be used for the conversion of the interferograms from geographical to disloc coordinates (km)
%

%  x_posting,y_posting is the distance in km between gridpoints measured in direction from the lower left corner (the origin ofthe  disloc system) 
%  For geographic coordinates we calculate an average x_posting,y_posting for the center of the interferogram. This approximation is neglible 
%  for small areas (e.g. volcanoes) but could be a problem for multi-frame interferograms at high latitudes. The variation of easting for 
%  1 degree difference at a given katitude is calcualted using the command                                 
%  lat=40;  (distance(lat,0,lat,1,almanac('earth','WGS84')) - distance(lat-1,0,lat-1,1,almanac('earth','WGS84'))) / distance(lat,0,lat,1,almanac('earth','WGS84')) * 100   
%  note that x_posting,y_posting is used for calculating coordinates after grid/quadtree deecomposition and for coordinates of the planes in the inversion   
%  It should be realtively easy to calculate a grid of proper coordinates (using distance) and use this through the inversion/modeling)                      
%  For UTM coordinates x_posting,y_posting follows directly from x_step,y_step (by dvision by 1000). For Geographic coordinates we calculate
%  FA Oct 2006 
%
if strcmp(basemap.x_unit,'degres') || strcmp(basemap.x_unit,'degrees')
    x_last    = basemap.x_first + basemap.x_step*size(basemap.data,2);
    y_last    = basemap.y_first + basemap.y_step*size(basemap.data,1);

    x_center  = (basemap.x_first+x_last)/2.;
    y_center  = (basemap.y_first+y_last)/2.;

    pt1       = [y_center,basemap.x_first] ;
    pt2       = [y_center,        x_last ] ;
    x_posting = distance(pt1,pt2,almanac('earth','WGS84')) / size(basemap.data,2) ;

    pt1       = [basemap.y_first,x_center] ;
    pt2       = [        y_last ,x_center] ;
    y_posting = distance(pt1,pt2,almanac('earth','WGS84'))/ size(basemap.data,1);
elseif strcmp(basemap.x_unit,'Meters')
    x_posting  = basemap.x_step/1000;          % convert into [km],
    y_posting  = abs(basemap.y_step)/1000;     % [km] and positive (because origin is now lower left),
else
    error('user error: x_unit %s not recognized -- exiting', basemap.x_unit) ;
end

basemap.x_posting=x_posting ;
basemap.y_posting=y_posting ;

%
% save data
%
  save(outfile,'basemap')
