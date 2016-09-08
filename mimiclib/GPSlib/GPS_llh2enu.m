function [GPSdist] = GPS_llh2enu(GPS);

%
% reads GPS structure and transform into distance changes 
%
% [GPSdist] = GPS_llh2enu(GPSlat);
%
% GPSlat is GPS structure from Read_GPS with coordinates change
%
% Use of GRS80 ellipsoide
%
% N. Gourmelen, February 2006
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GPSdist = GPS ;

for ni=1:length(GPS)

	latlength = distance(GPS(ni).lat(1),mean(GPS(ni).lon),GPS(ni).lat(size(GPS(ni).lat,1)),mean(GPS(ni).lon),almanac('earth','grs80')) ;
	lonlength = distance(mean(GPS(ni).lat),GPS(ni).lon(1),mean(GPS(ni).lat),GPS(ni).lon(size(GPS(ni).lon,1)),almanac('earth','grs80')) ;
    
	latratio  = latlength/abs(GPS(ni).lat(1)-GPS(ni).lat(size(GPS(ni).lat,1))) ; if isnan(latratio)  latratio = 0 ;  end
	lonratio  = lonlength/abs(GPS(ni).lon(1)-GPS(ni).lon(size(GPS(ni).lon,1))) ; if isnan(lonratio)  lonratio = 0 ;  end

	GPSdist(ni).n       = ( (GPS(ni).lat - mean(GPS(ni).lat) ) * latratio)*1000 ;
	GPSdist(ni).e       = ( (GPS(ni).lon - mean(GPS(ni).lon) ) * lonratio)*1000 ;
	GPSdist(ni).u       = GPS(ni).height - mean(GPS(ni).height)                 ;
	GPSdist(ni).error_e = GPS(ni).error_lon * lonratio * 1000                   ;
	GPSdist(ni).error_n = GPS(ni).error_lat * latratio * 1000                   ;
    GPSdist(ni).error_u = GPSdist(ni).error_height * 1000                       ;
	GPSdist(ni).lat     = mean(GPS(ni).lat)                                     ;
    GPSdist(ni).lon     = mean(GPS(ni).lon)                                     ;
    GPSdist(ni).height  = mean(GPS(ni).height)                                  ;
	GPSdist(ni).Unit    = 'm'                                                   ; 
end

GPSdist = rmfield(GPSdist,'error_lat') ;  GPSdist = rmfield(GPSdist,'error_lon') ;
