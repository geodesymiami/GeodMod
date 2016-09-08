function GPS = GPS_ITRF2local(GPS,opt)

%
% Rotate GPS vector and timeseries using poles from Sella et al., 2002 and method described by Cox and hart 1986
%
%  [outGPS] = GPS_ITRF2local(GPS,opt)
%
%    GPS: structure from Read_GPS
%
%      opt can be:
%
%        'plate': Which tectonic plate are your GPS on. Ideally should determine it from GPS location
%
%                   NorthAmerica (default)
%                   Nazca
%                   Pacific
%        
%    
% N. Gourmelen Feb - 2008
%

%% Deal options

defaultopt=struct(                                                             ...
         'plate'                               ,          'NorthAmerica'        ,            ...
    	 'method'     	                       ,             'linear'           )             ;            
    
    
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],[]); end ;
[opt]=process_defaultoptions(opt,defaultopt);  display(opt)
f=fieldnames(opt) ;
for i=1:length(f)  eval([char(f{i}) '= opt.(f{i}) ;' ]) ;  end


%% Compute predicted velocity at GPS site. ITRF97 Pole determination
switch(plate)
    
    
    
    case('Nazca')
        
        AngVel = 0.647e-6 / 180 * pi ;  radPole = [44.45 -99.49] / 180 * pi ;

        for ni=1:length(GPS)

            radGPS = [GPS(ni).lat GPS(ni).lon] / 180 * pi ;

            % Radius of the earth at GPS location
            R  = sqrt(( ((6378137^2 * cos(radGPS(1)))^2 + (6356752^2 * sin(radGPS(1)))^2) / ((6378137   * cos(radGPS(1)))^2 + (6356752   * sin(radGPS(1)))^2) ));

            % Convert to cartesian using spherical geometry
            [PoleX,PoleY,PoleZ] = ell2xyz(radPole(1),radPole(2),0,1,0) ;  [GPSX,GPSY,GPSZ] = ell2xyz(radGPS(1),radGPS(2),0,1,0) ;
            Polevect = [PoleX,PoleY,PoleZ] ;  GPSvect = [GPSX,GPSY,GPSZ] ;

            % Predicted velocity at the GPS location in cartesian coordiantes
            Vcart    = AngVel * R * cross(Polevect,GPSvect) ;

            % Transformation from cartesian to local neu coordinates
            T = [-sin(radGPS(1))*cos(radGPS(2))  -sin(radGPS(1))*sin(radGPS(2))  cos(radGPS(1)) ; ...
                 -sin(radGPS(2))                  cos(radGPS(2))                      0         ; ...
                 -cos(radGPS(1))*cos(radGPS(2))  -cos(radGPS(1))*sin(radGPS(2)) -sin(radGPS(1)) ];

            Vneu(:,ni) = T * Vcart' ;
            
        end
    
    case('NorthAmerica')
        
        AngVel = 0.199e-6 / 180 * pi ;  radPole = [-2.39 -79.08] / 180 * pi ;

        for ni=1:length(GPS)

            radGPS = [GPS(ni).lat GPS(ni).lon] / 180 * pi ;

            % Radius of the earth at GPS location
            R  = sqrt(( ((6378137^2 * cos(radGPS(1)))^2 + (6356752^2 * sin(radGPS(1)))^2) / ((6378137   * cos(radGPS(1)))^2 + (6356752   * sin(radGPS(1)))^2) ));

            % Convert to cartesian using spherical geometry
            [PoleX,PoleY,PoleZ] = ell2xyz(radPole(1),radPole(2),0,1,0) ;  [GPSX,GPSY,GPSZ] = ell2xyz(radGPS(1),radGPS(2),0,1,0) ;
            Polevect = [PoleX,PoleY,PoleZ] ;  GPSvect = [GPSX,GPSY,GPSZ] ;

            % Predicted velocity at the GPS location in cartesian coordiantes
            Vcart    = AngVel * R * cross(Polevect,GPSvect) ;

            % Transformation from cartesian to local neu coordinates
            T = [-sin(radGPS(1))*cos(radGPS(2))  -sin(radGPS(1))*sin(radGPS(2))  cos(radGPS(1))   ; ...
                 -sin(radGPS(2))                  cos(radGPS(2))                      0           ; ...
                 -cos(radGPS(1))*cos(radGPS(2))  -cos(radGPS(1))*sin(radGPS(2)) -sin(radGPS(1)) ] ;

            Vneu(:,ni) = T * Vcart' ;
            
        end
        
    case('Pacific')
        
        AngVel = 0.665e-6 / 180 * pi ;  radPole = [-64.21 112.74] / 180 * pi ;
        
        for ni=1:length(GPS)
            
            radGPS = [GPS(ni).lat GPS(ni).lon] / 180 * pi ;
            
            % Radius of the earth at GPS location
            R  = sqrt(( ((6378137^2 * cos(radGPS(1)))^2 + (6356752^2 * sin(radGPS(1)))^2) / ...
                ((6378137   * cos(radGPS(1)))^2 + (6356752   * sin(radGPS(1)))^2) ));

            %Ellipsoid = almanac('earth','GRS80') ;
            %[PoleX,PoleY,PoleZ] = ell2xyz(PoleLat,PoleLon,0,Ellipsoid(1)*1000,Ellipsoid(2)) ;
            %[GPSX,GPSY,GPSZ]    = ell2xyz(GPS(ni).lat,GPS(ni).lon,mean(GPS(ni).height),Ellipsoid(1)*1000,Ellipsoid(2)) ;
            
            % Convert to cartesian using spherical geometry 
            Ellipsoid = almanac('earth','GRS80') ;
            [PoleX,PoleY,PoleZ] = ell2xyz(radPole(1),radPole(2),0,Ellipsoid(1)*1000,Ellipsoid(2)^2)            ;  
            [GPSX,GPSY,GPSZ]    = ell2xyz(radGPS(1),radGPS(2),GPS(ni).height,Ellipsoid(1)*1000,Ellipsoid(2)^2) ;
            Polevect = [PoleX,PoleY,PoleZ] / sqrt(sum([PoleX,PoleY,PoleZ].^2)) ;  
            GPSvect  = [GPSX,GPSY,GPSZ]    / sqrt(sum([GPSX,GPSY,GPSZ].^2)) ;

            % Predicted velocity at the GPS location in cartesian coordiantes
            Vcart    = AngVel * R * cross(Polevect,GPSvect) ;

            % Transformation from cartesian to local eun coordinates
            T = [-sin(radGPS(1))*cos(radGPS(2))  -sin(radGPS(1))*sin(radGPS(2))  cos(radGPS(1))   ; ...
                 -sin(radGPS(2))                  cos(radGPS(2))                      0           ; ...
                 -cos(radGPS(1))*cos(radGPS(2))  -cos(radGPS(1))*sin(radGPS(2)) -sin(radGPS(1)) ] ;
            
            Vneu(:,ni) = T * Vcart' ;

        end
        
end

%% remove plate motion from GPS time series and/or rate

for ni=1:length(GPS)
    
    GPS(ni).n = GPS(ni).n - (GPS(ni).date_years - GPS(ni).date_years(1)) * Vneu(1,ni);
    GPS(ni).e = GPS(ni).e - (GPS(ni).date_years - GPS(ni).date_years(1)) * Vneu(2,ni);
    GPS(ni).u = GPS(ni).u - (GPS(ni).date_years - GPS(ni).date_years(1)) * Vneu(3,ni);
    
    GPS(ni).n = GPS(ni).n - mean(GPS(ni).n) ;
    GPS(ni).e = GPS(ni).e - mean(GPS(ni).e) ;
    GPS(ni).u = GPS(ni).u - mean(GPS(ni).u) ;
    
end



