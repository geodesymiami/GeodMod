function [InSARout,opt]=MakeProfile(data,opt);

%
% Project data matrix onto a profile. Output structure with two fields:
%
% InSARout.profileall has three vectors: [distance along profile from loca / value / dem (if present)]
% InSARout.profileav  has same vectors with averaged value if average is wanted
% InSARout.coordinates  is a 2*n vector with lat and long of points along profile
%
% Plot LOS in function of distance along the profile if plotres option is set to 'on'.
%
%
% [InSARout]=ProjectInSAR(data,opt)
%
%  opt=struct('datafield','datafit','projectpara',struct('profbounds',[30 30],'loca',[36.58 -117.67],'distmax',20,'dir',20));
% 
%
%   data      : igram structure
%
%
%   opt is
%
% 'projectpara'	 : Structure with parameters of the projection
%
%       loca        : Location of origine (fault?) [lat lon]
%       dir         : Direction in degree from north
%       profbounds  : Max length of profile from each side of loca.
%                     [North/East South/West]. North if 0<dir<45 or 135<dir<190; East otherwise
%                     Unit is km if proj=geo, pixel otherwise
%       distmax     : Maximum distance of the GPS from the profile 
%                     Unit is km if proj=geo, pixel otherwise
%       average     : Averaged the profile for smoothness (value in kilometers) (%FA 02/08 not working)
%       GridOneD    : number of datapoints in profile direction if given (similar to GridRows in SampelSar)
%                     the mean of all points within a distance bin is evaluated.
%
% 'dem'             : Structure with a dem for correlation
% 'dataformat'      : 'grid' (default) or 'vector' or 'gps'
% 'gpsfield'        : 'los_rate' (default), 'u_rate', 'e_rate', 'n_rate', 'horlos_rate'
% 'plotres'         : Plot if keyword present 
%
% N. Gourmelen, April 2007 
% FA, Feb 2008
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultopt=struct(                          ...
    'dem'                ,  'off'       ,   ...
    'dataformat'	     ,  'grid'      ,   ...
    'gpsfield'           ,  'los_rate'  ,   ...
    'datafield'          ,  'data'      ,   ...
    'GridOneD'           ,  'off'       ,   ...
    'plotres'            ,  'off'   );
defaultopt.projectpara=struct(       ...
    'DoIt'               ,  'on'    ,   ...
    'GridOneD'           ,  'off'   );

if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],[]); end ;
[opt]=process_defaultoptions(opt,defaultopt);  
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get distance along profile and remove value outside of range distance from the profile %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Convert DEM to Vector
if isstruct(dem)  DEMvect=dem.data(:);  end

%Check if Grid or Vector Data

if strcmp(dataformat, 'grid') == 1
    % Vectorise InSAR data
    InSARvect = data.(genvarname(datafield))(:) ;
    % Build coordinates matrice
    outind = [] ;

    lonvect = [data.x_first:data.x_step:data.x_first+data.x_step*size(data.data,2)]  ;  %FA 3/14 changed from width, filelength to size, may have to switch index
    latvect = [data.y_first:data.y_step:data.y_first+data.y_step*size(data.data,1)] ;
    
    [InSARlon,InSARlat] = meshgrid(lonvect,latvect) ;
    
elseif strcmp(dataformat, 'vector') == 1
    InSARvect    = data.(genvarname(datafield)) ;
    InSARvectori = InSARvect                    ;

    % Build coordinates matrice
    outind = [] ;

    InSARlon = data.lon ;  InSARlonori = InSARlon ;
    InSARlat = data.lat ;  InSARlatori = InSARlat ;

elseif strcmp(dataformat, 'gps') == 1
    
    InSARlon = [data.lon] ;  InSARlat = [data.lat] ;
    
    if isfield(data(1),gpsfield)
        InSARvect = [data.(genvarname(gpsfield))] ;                                                                      ;
    end
    
end


% Shift origin to the input location

  InSARlon_shift=InSARlon(:)-projectpara.loca(2);InSARlat_shift=InSARlat(:)-projectpara.loca(1);

% Rotate coordinates by -projectpara.dir

  InSARlonrot=InSARlon_shift*cos(-projectpara.dir/180*pi)+InSARlat_shift*sin(-projectpara.dir/180*pi);
  InSARlatrot=-InSARlon_shift*sin(-projectpara.dir/180*pi)+InSARlat_shift*cos(-projectpara.dir/180*pi);

% Shift back origin so coordinates

  newlon=projectpara.loca(2)+InSARlonrot;newlat=projectpara.loca(1)+InSARlatrot;

% Project on profile - To Get Distance From Profile 

  newlonshift=0;

% Shift in pseudo-latitude

  newlatshift=newlat-projectpara.loca(1);

% Rotate by projectpara.dir

  InSARlonrotonprof=newlonshift*cos(projectpara.dir/180*pi)+newlatshift*sin(projectpara.dir/180*pi);
  InSARlatrotonprof=-newlonshift*sin(projectpara.dir/180*pi)+newlatshift*cos(projectpara.dir/180*pi);

% Shift back
  
  newlononprofile=InSARlonrotonprof+projectpara.loca(2);newlatonprofile=InSARlatrotonprof+projectpara.loca(1);

% Distance of original coordinates from projected ones on profile 

data.proj='geo';      % FA 3/2014  Not sure why Noel uses geo and what it means
 if strcmp(data.proj,'geo') 
     distfromprofile = distance(newlatonprofile,newlononprofile,InSARlat(:),InSARlon(:),almanac('earth','grs80')) * 1000;
 else
     distfromprofile = sqrt( (newlatonprofile-InSARlat(:)).^2 + (newlononprofile-InSARlon(:)).^2 ) ;
 end
% Find and remove points farther than projectpara.distmax

  if isfield(projectpara,'distmax')
      notok=find(distfromprofile > projectpara.distmax);
  end

% Distance from projectpara.loca on profile

data.proj='geo';   %FA 3/14  Not sure why Noel uses "geo"
if strcmp(data.proj,'geo')  
      distonprofile = distance(newlatonprofile,newlononprofile,projectpara.loca(1),projectpara.loca(2),almanac('earth','grs80')) * 1000;
  else
      distonprofile = sqrt( (newlatonprofile-projectpara.loca(1)).^2 + (newlononprofile-projectpara.loca(2)).^2 ) ;
  end
  
% Get sign of Distance. West and South negative

  if (0<=projectpara.dir & projectpara.dir<=45) | (135<=projectpara.dir & projectpara.dir<=190)
      signind = find(newlatonprofile<projectpara.loca(1)) ;
      distonprofile(signind) = -distonprofile(signind)    ;
  else
      signind = find(newlononprofile<projectpara.loca(2)) ;
      distonprofile(signind) = -distonprofile(signind)    ;
  end

% Find and remove points farther than projectpara.profbounds

if isfield(projectpara,'profbounds')
    if length(projectpara.profbounds) ~= 2  projectpara.profbounds(2) = projectpara.profbounds(1) ; end
    
    upno   = find(distonprofile > projectpara.profbounds(1))  ;
    downno = find(distonprofile < -projectpara.profbounds(2)) ;

    % Find NaNs

    nn = find(isnan(InSARvect)==1) ;

    % Remove data-points
%FA 3/2014. With the parameter I have it removed everything. This needs t
%be checked
    toremove                  = unique([nn;notok;upno;downno]) ;
    newlononprofile(toremove) = []                  ;
    newlatonprofile(toremove) = []                  ;
    distonprofile(toremove)   = []                  ;
    InSARvect(toremove)       = []                  ;

end
  
% Sort by distance from the fault

  [distonprofile,i_distonprofile]  = sort(distonprofile)              ;
  newlononprofile                  = newlononprofile(i_distonprofile) ;
  newlatonprofile                  = newlatonprofile(i_distonprofile) ;
  InSARvect                        = InSARvect(i_distonprofile)       ;
  
  if isstruct(dem)  DEMvect(toremove) = [] ; DEMvect = DEMvect(i_distonprofile);  end

  if isfield(projectpara,'average')
      avInSARvect = smooth(InSARvect,distonprofile,projectpara.average);
  else
      avInSARvect = InSARvect;  %avdistonprofile=distonprofile;
  end

% Prepare profile area

  InSARlon(toremove) = [] ;  InSARlat(toremove) = [] ;  dd = unique(InSARlat) ;
  InSARlon = InSARlon(i_distonprofile) ;  InSARlat = InSARlat(i_distonprofile) ;

  InSARlonrotonprof = newlonshift*cos(projectpara.dir/180*pi)+newlatshift*sin(projectpara.dir/180*pi)  ;
  InSARlatrotonprof = -newlonshift*sin(projectpara.dir/180*pi)+newlatshift*cos(projectpara.dir/180*pi) ;

% Build box

  unitlon      = sin(projectpara.dir/180*pi) ;
  unitlat      = cos(projectpara.dir/180*pi) ;
  unitlonortho = -unitlat                    ;
  unitlatortho = -unitlon                    ;

  if strcmp(data.proj,'geo') 
      AlongTrackunit  = distance(projectpara.loca(1),projectpara.loca(2),projectpara.loca(1)+unitlat,projectpara.loca(2)+unitlon,almanac('earth','grs80'))           * 1000 ;
      AcrossTrackunit = distance(projectpara.loca(1),projectpara.loca(2),projectpara.loca(1)+unitlatortho,projectpara.loca(2)+unitlonortho,almanac('earth','grs80')) * 1000 ;
  else
      AlongTrackunit  = sqrt( (unitlat).^2 + (unitlon).^2 ) ;
      AcrossTrackunit = sqrt( (unitlatortho).^2 + (unitlonortho).^2 ) ;
  end
  
  AlongTrackFactN = projectpara.profbounds(1)/AlongTrackunit ;
  AlongTrackFactS = projectpara.profbounds(2)/AlongTrackunit ;

  AcrossTrackFact = projectpara.distmax/AcrossTrackunit      ;

  lon_ul = -(AcrossTrackFact)*cos(projectpara.dir/180*pi)+(AlongTrackFactN)*sin(projectpara.dir/180*pi) ;
  lon_ur = (AcrossTrackFact)*cos(projectpara.dir/180*pi)+(AlongTrackFactN)*sin(projectpara.dir/180*pi)  ;
  lon_ll = -(AcrossTrackFact)*cos(projectpara.dir/180*pi)-(AlongTrackFactS)*sin(projectpara.dir/180*pi) ;
  lon_lr = (AcrossTrackFact)*cos(projectpara.dir/180*pi)-(AlongTrackFactS)*sin(projectpara.dir/180*pi)  ;

  lat_ul = (AcrossTrackFact)*sin(projectpara.dir/180*pi)+(AlongTrackFactN)*cos(projectpara.dir/180*pi)  ;
  lat_ur = -(AcrossTrackFact)*sin(projectpara.dir/180*pi)+(AlongTrackFactN)*cos(projectpara.dir/180*pi) ; 
  lat_ll = (AcrossTrackFact)*sin(projectpara.dir/180*pi)-(AlongTrackFactS)*cos(projectpara.dir/180*pi)  ; 
  lat_lr = -(AcrossTrackFact)*sin(projectpara.dir/180*pi)-(AlongTrackFactS)*cos(projectpara.dir/180*pi) ;

% Polygon

  opt.projectpara.polygon=[[lat_ul lat_ur lat_lr lat_ll lat_ul]'+projectpara.loca(1) [lon_ul lon_ur lon_lr lon_ll lon_ul]'+projectpara.loca(2)]

% Quick and dirty 1D gridding/averaging if projectpara.GridOneD given , FA Feb 2008

  GridOneD = projectpara.GridOneD;

  if GridOneD
    
     InSARvect = resample(double(InSARvect),1,round(length(distonprofile)/GridOneD),round(length(distonprofile)/GridOneD)) ;
     distonprofile = resample(double(distonprofile),1,round(length(distonprofile)/GridOneD),round(length(distonprofile)/GridOneD)) ;
     
  else
      InSARout.profileav    = [distonprofile(:) avInSARvect(:) ] ;   
      InSARout.coordinates  = [newlatonprofile newlononprofile]  ;
  end

% Assign output Structure

  InSARout.profileall   = [distonprofile(:) InSARvect(:)   ] ;


