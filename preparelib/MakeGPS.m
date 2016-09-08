function [GPSfile,GPS]=MakeGPS(makegpsopt)
%MakeGPS     - reads GPS data in ascii format and writes as *mat file. 
%
%  usage:  [GPSfile]    =MakeGPSData(makegpsopt)
%          [GPSfile,GPS]=MakeGPSData(makegpsopt)
%  
%  input:   makegpsopt constains the following fields:
%
%           GPStxtfile - name of ascii file
%           GPSformat  - foramt of ascii file (currently only supports 'BenBrooks','INGV-CT')
%           RemoveList - stations to exclude. Not yet implemented
%           GPSfile    - name of mat file
%  output:  GPSfile is the name of the saved *mat file containing the GPS structure
%
%  TODO: Needs to be properly updated to Geodmod format
%  TODO: We probably should have a substructure for the GPS data (e.g. makegpsopt.format, 
%        makegpsopt.infile, makegpsopt.removelist etc.)
%
%  Falk Amelung, June 2007
global dir_out
defaultopt=struct(                                                         ...
        'DoIt'               ,        'on'                    ,            ...
        'DoPlot'             ,        'on'                    ,            ...
        'GPStxtfile'         ,        'off'                   ,            ...
        'GPSfile'            ,        'off'                   ,            ...
        'filetype'           ,        'rates'                 ,            ...
        'GPSformat'          ,        'BenBrooks'             ,            ...
        'radarCoordinates'   ,        'off'                   ,            ...
        'GPSsta'             ,        'off'      );    
[makegpsopt]=process_defaultoptions(makegpsopt,defaultopt);  
f=fieldnames(makegpsopt) ; for i=1:length(f) eval([char(f{i}) '= makegpsopt.(f{i}) ;' ]) ; end
if  ~DoIt   return; end
if  ~GPStxtfile(1) && ~GPSfile(1) GPS=0;   return; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if GPStxtfile(1) CheckInOut(GPStxtfile,''); end             % check whether GPStxtfile exist
if ~GPSfile(1)                                              % name the GPS*mat file
    [pathstr,name,ext] = fileparts(GPStxtfile);
    GPSfile            = [fullfile(dir_out,name) '.mat'];
end
if  CheckInOut('',GPSfile)  load(GPSfile); return; end      % check whether GPStxtfile exist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% NOEL %%%%
fid = fopen(GPStxtfile);
%
%if ~any(strmatch(GPSformat,{'BenBrooks' 'INGV-CT' 'UNR-NOUP'}))
%   errordlg(sprintf('GPSformat not recognized %s',GPSformat)) ;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%opt = struct('source',GPSformat); 
%GPS = Read_GPS(GPStxtfile,opt);

% Read GPS
%opt = struct('GPSformat',GPSformat,'coordSystem',coordSystem,'plate',plate,'filetype',filetype);

if ~isfield(makegpsopt,'source') makegpsopt.source=makegpsopt.GPSformat; end  % this is necessary because sometimes source and sometimes GPSformat is used. FA 2/2010

GPS = Read_GPS(GPStxtfile,makegpsopt) ;

%% Prepare data when time series
if sum(strcmp(filetype,{'tmsr','tmsrBaseline'}))
    if strcmp(coordSystem,'ITRF') 
        for ni=1:length(GPS)
            GPS(ni).X = GPS(ni).lon    ;
            GPS(ni).Y = GPS(ni).lat    ;
            GPS(ni).Z = GPS(ni).height ;
            GPS(ni).error_X = GPS(ni).error_lon    ;
            GPS(ni).error_Y = GPS(ni).error_lat    ;
            GPS(ni).error_Z = GPS(ni).error_height ;
            
            [lat,lon,alt] = ecef2lla(GPS(ni).X,GPS(ni).Y,GPS(ni).Z) ;
            [e_lat,e_lon,e_alt] = ecef2lla(GPS(ni).X+GPS(ni).error_X,GPS(ni).Y+GPS(ni).error_Y,GPS(ni).Z+GPS(ni).error_Z) ;
            
            GPS(ni).lon = rad2deg(lon) -360 ;
            GPS(ni).lat = rad2deg(lat) ;
            GPS(ni).height = alt       ;
            GPS(ni).error_lon = abs(rad2deg(lon) - rad2deg(e_lon)) ;
            GPS(ni).error_lat = abs(rad2deg(lat) - rad2deg(e_lat)) ;
            GPS(ni).error_height = abs(alt - e_alt) ;
            
        end
    end
    
    GPS = GPS_llh2enu(GPS)                                          ;
    %GPS = GPS2local(GPS,struct('plate',inputdataopt.gpsplate)) ;
    
    if isfield(makegpsopt,'radarlook')  GPS = gps2los(GPS,makegpsopt.radarlook) ;  end
    GPS = GPS_tmsr2rate(GPS)    ;
                  
    %GPS = Geo2Radar(GPS,struct(inputdataopt))                       ;
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%gps_cov=eye(length(sites)*3) ;
%gps_sig=gps_sig/1000;
% May 1 2007: Something is not right with the covarianec/sigma. Need to figure this out.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Remove stations out of area %%%%%%%%%%%%%%%%%%%%%%
removelist={'MLSP' '75FL' 'VO92' 'MOKP' 'V134' '133T'};
removelist=[];  torem=[];
for ni=1:length(removelist)
    for ny=1:length(GPS)
        if strmatch(GPS(ny).station,removelist(ni))  torem=[torem;ny];  end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Remove stations out of DEM area %%%%%%%%%%%%%%%%%%
y_first  = basemap.y_first;
x_first  = basemap.x_first;
y_last   = basemap.y_first + basemap.y_step*size(basemap.data,1);
x_last   = basemap.x_first + basemap.x_step*size(basemap.data,2);
ind      = find([GPS.lat]>y_first | [GPS.lat]<y_last);
GPS(ind) = [];                                       
ind      = find([GPS.lon]>x_last  | [GPS.lon]<x_first);
GPS(ind) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Transform coordinates of gps stations into UTM and more %%%
if isstruct(radarCoordinates)
    inStations = getfields(GPS,'station','cell') ;  rCoordstations = fieldnames(radarCoordinates) ;
    for ny=1:length(rCoordstations)
        dI(ny) = find(strcmp(inStations,rCoordstations(ny))==1) ;
        yx     = radarCoordinates.(rCoordstations{ny})          ;
        GPS(dI(ny)).y = yx(1) ;  GPS(dI(ny)).x = yx(2) ;  GPStmp(ny) = GPS(dI(ny)) ;
    end
    GPS = GPStmp ;
else
    if ~isempty(GPS)                                         % FA 1/2011 Added so that it does not fail when no GPS in area (e.g. Wells)
        zone = utmzone( GPS(1).lat, GPS(2).lon );
       figure;axesm('mapprojection','utm','zone',zone);      % Here the Fig Axes is only opened for the map projection
    end
    for ni=1:length(GPS)
        [GPS(ni).utm_x,GPS(ni).utm_y] = mfwdtran( GPS(ni).lat, GPS(ni).lon  );
        GPS(ni).Unit                  = 'm/yr';
        GPS(ni).station               = char(GPS(ni).station);
        xy                            = lola2xy([GPS(ni).lon;GPS(ni).lat],basemap,1);
        GPS(ni).x                     = xy(1);
        GPS(ni).y                     = xy(2);
    end % UTM Easting and Northing
    close
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(GPSfile,'GPS')    % same format as SJonni for hector, I hope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (DoPlot && ~isempty(GPS))
    lats     = getfields(GPS,'lat');
    lons     = getfields(GPS,'lon');
    stations = getfields(GPS,'station','cell');
    east     = getfields(GPS,'e_rate');
    north    = getfields(GPS,'n_rate');
    up       = getfields(GPS,'u_rate');
    utm_x    = getfields(GPS,'utm_x');
    utm_y    = getfields(GPS,'utm_y');
    enu      = [east';north';up'];
    
    latlim = [basemap.y_first+basemap.y_step*basemap.file_length,basemap.y_first];
    lonlim = [basemap.x_first,basemap.x_first+basemap.x_step*basemap.width];

    %%% Plot GPS Data, Yunjun, 2015-12-02
    scsz = get(0,'ScreenSize');
    hscale = max(sqrt(east.^2+north.^2)); vscale = max(up);
    figure('Position',[10 10 0.5*scsz(4) 0.8*scsz(4)],'Name','GPS Data');

    h1=subplot(2,1,1); hold on;                                                % plot Horizontal deformation
    image(lonlim,latlim,basemap.shade); hold on;
    plot(lons,lats,'s','MarkerSize',6);
    text(lons*1.00004,lats,stations,'FontSize',6);
    quiver(lons,lats,east/hscale,north/hscale,'b','LineWidth',2);
    title('Horizontal deformation'); box on;
    legend('GPS sites','Data'); legend('boxoff');
    axis equal; axis tight;
    
    h2=subplot(2,1,2); hold on;                                                % plot Vertical deformation
    image(lonlim,latlim,basemap.shade); hold on;
    plot(lons,lats,'s','MarkerSize',4);
    text(lons*1.00004,lats,stations,'FontSize',6);
    quiver(lons,lats,zeros(size(up)),up/vscale,'b','LineWidth',2);
    title('Vertical deformation'); box on;
    axis equal; axis tight;
    legend('GPS sites','Data'); legend('boxoff')
    % link lim for x/y
    xlim = [min([get(h1,'Xlim'),get(h2,'Xlim')]),max([get(h1,'Xlim'),get(h2,'Xlim')])]; 
    ylim = [min([get(h1,'Ylim'),get(h2,'Ylim')]),max([get(h1,'Ylim'),get(h2,'Ylim')])]; 
    xwid = xlim(2)-xlim(1);             xlim = xlim+[xwid*-0.05 xwid*0.05];
    ywid = ylim(2)-ylim(1);             ylim = ylim+[ywid*-0.05 ywid*0.05];
    set(h1,'Xlim',xlim,'Ylim',ylim);    set(h2,'Xlim',xlim,'Ylim',ylim);
    % save fig to pdf file
    if  ~CheckInOut('',[ dir_out filesep 'plot' filesep 'GPSData.pdf' ])
        logplot('',[dir_out filesep 'GPSData'],'');
    end

%     figure
%     plot  (lons,lats,'d');
%     text  (lons,lats,stations,'FontSize',6,'VerticalAlignment','bottom','Color', 'k' );
%     [H,H_title]=DM_Quiver([ lons';lats' ],enu(:),GPS(1).cov,0.00001*1000) ; set(H,'color','b')
%     title(H_title);
%     %figure
%     %plot  (utm_x,utm_y,'d');
%     %text  (utm_x,utm_y,stations,'FontSize',6,'VerticalAlignment','bottom','Color', 'k' );
%     %H=DM_Quiver([ utm_x';utm_y' ],enu(:),GPS(1).cov,1.0*1000) ; set(H,'color','b')

end
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if GPStxtfile(1) || GPSfile(1)
%   [GPS] = MakeGPSData(opt) ;
%   if strcmp(basemap.x_unit,'Meters')         GPSsta = [GPS.utm; GPS.enu] ;            
%      elseif strcmp(basemap.x_unit,'degres')  GPSsta = [GPS.lola;  GPS.enu] ;
%      else error ('x_unit not recognized: %s -- exiting',basemap.x_unit)      
%   end
%end
