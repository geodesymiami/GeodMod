function [satparameter]=extract_hardwired_satparameters(in_name, parameter)
%extract_hardwired_satparameters  - returns satellite parameter (e.g. LOSvector, wavelength ) 
%
%usage: [satparameters]=extract_hardwired_satparameters(in_name,parameter);
%
%Input:   in_name       file or path name containing satellite name and beam (e.g. RsatD1, EnvD2) 
%         parameter     parameter to extract (currently works for 'sat_height','wavelength','LOSvector')
%
%Output:  satparameter  value for 'wavelength','sat_height' or 'LOSvector'
%
%  uses the extract_ProjectName function to extract the satellite from the filename
%
%  ATTENTION: azimuth is very inaccurate, could be improved
%
%  Falk Amelung, September 2006
%  FA Oct 5 2006  inverted sign of vertical component so that it is positive
%  FA Jan 2008    renamed from GenerateLOSvector. Introduced extract_Projectname and 'wavelength' and 'sat_height' functionality.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isstruct(in_name)  if isfield(in_name,'heading') azimuth = in_name(1).heading/180*pi ;  end ;  in_name = in_name(1).sat ;  end

[junk,junk,satbeam,sat]= extract_ProjectName(in_name);

switch parameter
    case {'wavelength'}
        switch sat
            case{'Rsat','Ers','Env'};  wavelength = 0.056 ;
            case{'Alos','Jers'};       wavelength = 0.236 ;
            case{'Csk'};       wavelength = 0.0312283810416667 ;  %Anieri 5/16
        end
        satparameter=wavelength;
    case {'sat_height'}
        switch sat
            case{'Ers','Env'};        sat_height = 785000 ;
            case{'Rsat'};             sat_height = 793000 ;
            case{'Jers'};             sat_height = 574000 ;
            case{'Alos'};             sat_height = 700000 ;
            case{'Csk'};             sat_height = 625844.2166 ;  %Anieri 5/16
        end
        satparameter=sat_height;
    case {'LOSvector','incangle','azimuth'}
        switch satbeam
            case {'RsatD1','RsatSD1'};  incangle = 23.5/180*pi ; if ~exist ('azimuth','var') azimuth =  -169/180*pi ;  end
            case {'RsatD2','RsatSD2'};  incangle = 27.7/180*pi ; if ~exist ('azimuth','var') azimuth =  -169/180*pi ;  end
            case {'RsatD4','RsatSD4'};  incangle = 37.0/180*pi ; if ~exist ('azimuth','var') azimuth =  -169/180*pi ;  end
            case {'RsatD6','RsatSD6'};  incangle = 43.5/180*pi ; if ~exist ('azimuth','var') azimuth =  -169/180*pi ;  end
            case {'RsatD7','RsatSD7'};  incangle = 47.2/180*pi ; if ~exist ('azimuth','var') azimuth =  -169/180*pi ;  end
            case {'RsatA3','RsatSA3'};  incangle = 30.5/180*pi ; if ~exist ('azimuth','var') azimuth =   -11/180*pi ;  end
            case {'RsatA5','RsatSA5'};  incangle = 39.4/180*pi ; if ~exist ('azimuth','var') azimuth =   -11/180*pi ;  end
            case {'RsatA6','RsatSA6'};  incangle = 43.5/180*pi ; if ~exist ('azimuth','var') azimuth =   -11/180*pi ;  end
            case {'RsatA7','RsatSA7'};  incangle = 47.2/180*pi ; if ~exist ('azimuth','var') azimuth =   -11/180*pi ;  end
                
            case {'ErsD'};              incangle = 23.0/180*pi ; if ~exist ('azimuth','var') azimuth =  -169/180*pi ;  end
            case {'ErsA'};              incangle = 23.0/180*pi ; if ~exist ('azimuth','var') azimuth =   -11/180*pi ;  end
            case {'ErsDmai'};           incangle = 90.0/180*pi ; if ~exist ('azimuth','var') azimuth =   101/180*pi ;  else azimuth = mod(azimuth,pi)+pi/2 ;  end 
            case {'ErsAmai'};           incangle = 90.0/180*pi ; if ~exist ('azimuth','var') azimuth =  -101/180*pi ;  else azimuth = azimuth - pi/2        ;  end 
                
            case{'EnvD2'};              incangle = 23.0/180*pi ; if ~exist ('azimuth','var') azimuth =  -169/180*pi ;  end
            case{'EnvD6'};              incangle = 41.0/180*pi ; if ~exist ('azimuth','var') azimuth =  -169/180*pi ;  end
            case{'EnvA2'};              incangle = 23.0/180*pi ; if ~exist ('azimuth','var') azimuth =   -11/180*pi ;  end
            case{'EnvA6'};              incangle = 41.0/180*pi ; if ~exist ('azimuth','var') azimuth =   -11/180*pi ;  end
            case{'EnvD7'};              incangle = 43.9/180*pi ; if ~exist ('azimuth','var') azimuth =  -169/180*pi ;  end   

            case{'CskA'};              offnadir = 33.6/180*pi ; if ~exist ('azimuth','var') azimuth =  -11.7/180*pi ;  end   %Anieri 5/16 verify offnadir = average lookangle
            case{'CskD'};              offnadir = 34.8/180*pi ; if ~exist ('azimuth','var') azimuth = -168.3/180*pi ;  end   %Anieri 5/16 verify offnadir = average lookangle
                
            case{'JersD'};              offnadir = 35.0/180*pi ; if ~exist ('azimuth','var') azimuth =  -169/180*pi ;  end
            case{'AlosD'};              offnadir = 34.3/180*pi ; if ~exist ('azimuth','var') azimuth =-167.8/180*pi ;  end
            case{'AlosD04'};            offnadir = 25.8/180*pi ; if ~exist ('azimuth','var') azimuth =-167.8/180*pi ;  end
            case{'AlosA'};              offnadir = 34.3/180*pi ; if ~exist ('azimuth','var') azimuth = -12.2/180*pi ;  end
            otherwise
                errordlg(sprintf('satbeam not recognized %s',satbeam)) ;
        end
        
        % convert from look angle to incidence angle. %jb nov07
        
        if ~exist('incangle','var') 
              rearth       = 6.36e6;                                                                   % hard-wired parameter in m
              sat_height   = extract_hardwired_satparameters(in_name,'sat_height');
              incangle     = asin((sat_height+rearth)/rearth *sin(offnadir)) ;
         end
        radarlook    = [-sin(incangle)*cos(azimuth) sin(incangle)*sin(azimuth) cos(incangle)];   % Nov 25 2009 change
        
        switch parameter
            case {'LOSvector'};   satparameter = radarlook; 
            case {'incangle'};    satparameter = rad2deg(incangle);  
            case {'azimuth'};     satparameter = rad2deg(azimuth) + 90;   
        end
    otherwise
        errordlg(sprintf('parameter not recognized %s',parameter)) ;
        end
        %     FA 1/2011" ATTENTION: The azimuth is probably wrong. It is different from the heading in roi_pac's geo_incidence.unw
        %     FA 1/2011: There used to be a confusion. For all satellites
        %     the lookangle was given but it was infact the incidence
        %     angle, except for Alos in which case the off-nadir angle (or
        %     look angle was given)
        %    ALOS:  off-nadir=34.3:    http://earth.esa.int/download/alos/PALSAR_info_users_v1.1.pdf
        %    RSAT-1:                   http://gs.mdacorporation.com/partners/software/d4_3-0-1.zip    
        %                        S1    S2    S3    S4    S5    S6    S7
        %    incidence angle  20-27 24-31 30-37 34-40 36-42 41-46 45-49
        %    mean incangle     23.5  27.5  33.5  37.0  39.0  43.5  47.0
        %                             
        %                       F1     F2    F3    F4    F5
        %                     37-40 39-42 41-44 43-46 45-48
        %                      38.5  40.5  42.5  44.5  46.5
        
