function [GPS] = Read_GPS(datafiles,opt);

%
% reads all the GPS and saves all the GPS measurements relevant to the study area in a file called GPS.mat
%
% [GPS] = Read_GPS(datafile,struct('filetype','tmsr','source','Miami','Unit','mm'));
%
%    
%   datafiles     :   path relevant to the files containing the GPS measurement
%                   if rates, the first line should contain keywords "lat lon
%                   vn ve vu sn se su corr sta" in the order of the info in the
%                   file.
%                   By default, it reads the order above.
%
% opt is:
%
% filetype      :  'rates' or 'tmsr'. Type of file in input, time series or rates (default).
%
% parameters    :  Vector with keywords: [lat lon vn ve vu sn se su corr]
%
% Unit	        :  'm ; cm ; mm or llm (lat-lon-meter)'. assumes 'm'. Carefull that in llm, the height has to be in meter!!!
%
% source	:   'SOPAC','Miami', 'BenBrooks', 'INGV-CT', ...
%
%
% BoundingBox   :   [x_min x_max y_min y_max] are the geographic coordinates of the study area
%
% plotgps whether or not you want the timeseries to be plotted
%
%
% M. Manzo & N. Gourmelen, February 2006
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultopt=struct(                                                  ...
    'filetype'                  ,        'rates'              ,       ...
    'Unit'			            ,	     'off'		        ,   	...
    'source'			        ,	     'default'		        ,	    ...
    'BoundingBox'               ,        'off'              ,       ...
    'parameters'                ,        'off'              ,       ...
    'plotgps'                   ,        'off'      )                ;


%if strcmp(datafiles(end-2:end),
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],[]); end ;
[opt]=process_defaultoptions(opt,defaultopt);  %display(opt)
f=fieldnames(opt) ;
for i=1:length(f)
    eval([char(f{i}) '= opt.(f{i}) ;' ]) ;
end

%%%%%%%%%%%%%%%%%%%
% Check variables %
%%%%%%%%%%%%%%%%%%%

%if ~filetype  filetype='rates';  end
if strcmp(datafiles(end-3:end),'.mat')  filetype='matfile';  end
%if ~source    source='default';  end

if ~Unit
    unitsf=1;
elseif strcmp(Unit,'llm')
    unitsf=1;
elseif strcmp(Unit,'m')
    unitsf=1;
elseif strcmp(Unit,'cm')
    unitsf=0.01;
elseif strcmp(Unit,'mm')
    unitsf=0.001;
else
end

%%%%%%%%%%%%%%%%%%%
% Reads the files %
%%%%%%%%%%%%%%%%%%%

filetype
source

switch(filetype)

    case('rates')

        switch(source)

	    case('FinnurPalsson')
                fid = fopen(datafiles,'r');
                C   = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f%f','headerLines',1,'Delimiter',',') ;
                fclose(fid) ;
                datafile2 = input('Name of datafile: ','s') ;
                fid = fopen(datafile2) ;
                C2  = textscan(fid,'%s%f%f%f%f%f%f%f%f%f','headerLines',1) ;
                fclose(fid) ; 
                stations = C{1} ;  lat = C{6}+C{7}/60 ;  lon = C{8}+C{9}/60 ;  height = C{10} ;  
                vel = C2{10} ;  heading = C2{8} ;  
                for ni=1:size(stations,1)
                    ny = find(strcmp(C2{1},stations(ni))==1) ;
                    if size(vel(ny),1)==2 ;  vel(ny)  = mean(vel(ny)) ;  heading(ny) = mean(heading(ny)) ;  ny = ny(1) ;  end
                    GPS(ni)=struct('station',stations(ni),'lat',lat(ni),'lon',lon(ni),'vel',vel(ny),'heading',heading(ny),'height',height(ni));
                end
                
            case('SOPAC')

            case('Miami')

            case('BenBrooks')
                fid=fopen(datafiles,'r');
                C=textscan(fid,'%s%f%f%f%f%f%f%f%f','headerLines',1) ;
                stations=C{1};  lat=C{2}; lon=C{3};  east=C{6};  north=C{4};  up=C{8};  error_east=C{7};  error_north=C{5};  error_up=C{9};
                for ni=1:size(stations,1)
%                     GPS(ni)=struct('station',stations(ni),'lat',lat(ni),'lon',lon(ni),'e_rate',east(ni)/1000,'n_rate',north(ni)/1000, ...
%                         'u_rate',up(ni)/1000,'e_error',error_east(ni)/1000,'n_error',error_north(ni)/1000,'u_error',error_up(ni)/1000 ...
%                         ,'cov',zeros(size(stations,1)*3,size(stations,1)*3));
                   GPS(ni)=struct('station',stations(ni),'lat',lat(ni),'lon',lon(ni),'e_rate',east(ni).*unitsf,'n_rate',north(ni).*unitsf, ...
                        'u_rate',up(ni).*unitsf,'e_error',error_east(ni).*unitsf,'n_error',error_north(ni).*unitsf,'u_error',error_up(ni).*unitsf ...
                        ,'cov',zeros(size(stations,1)*3,size(stations,1)*3)); 
                    
                end
            case('INGV-CT')
                fid=fopen(datafiles,'r');
                C=textscan(fid,'%s%f%f%f%f%f%f%f%f%f','headerLines',1) ;
                stations=C{1};  utm_y=C{3}; utm_x=C{2};  east=C{5};  north=C{6};  up=C{7};  error_east=C{8};  error_north=C{9};  error_up=C{10};
                zone = '33N'; axesm('mapprojection','utm','zone',zone);  [lat,lon]=minvtran( utm_x, utm_y );   % UTM Easting and Northing
                for ni=1:size(stations,1)
                    GPS(ni)=struct('station',stations(ni),'lat',lat(ni),'lon',lon(ni),'utm_x',utm_x(ni),'utm_y',utm_y(ni),'e_rate',east(ni)/1000, ...
                        'n_rate',north(ni)/1000,'u_rate',up(ni)/1000,'e_error',error_east(ni)/1000,'n_error',error_north(ni)/1000, ...
                        'u_error',error_up(ni)/1000,'cov',zeros(size(stations,1)*3,size(stations,1)*3));
                end

            case('UNR-NOUP')
                fid=fopen(datafiles,'r');
                C=textscan(fid,'%f%f%f%f%f%f%s','headerLines',0) ;
                stations=C{7};  lat=C{2}; lon=C{1};  east=C{3};  north=C{4};  up=C{1}*0;  error_east=C{5};  error_north=C{6};  error_up=C{1}*0;
                for ni=1:size(stations,1)
                    GPS(ni)=struct('station',stations(ni),'lat',lat(ni),'lon',lon(ni),'e_rate',east(ni)/1000,'n_rate',north(ni)/1000, ...
                        'u_rate',up(ni)/1000,'e_error',error_east(ni)/1000,'n_error',error_north(ni)/1000,'u_error',error_up(ni)/1000 ...
                        ,'cov',zeros(size(stations,1)*3,size(stations,1)*3));
                end
            case('EricCalais1')
                fid=fopen(datafiles,'r');
                C=textscan(fid,'%f%f%f%f%f%f%f%s%f%f','headerLines',1) ;
                stations=C{8};  lat=C{2}; lon=C{1};  east=C{3};  north=C{4};  up=C{9};  error_east=C{5};  error_north=C{6};  error_up=C{10};
                for ni=1:size(stations,1)
                    GPS(ni)=struct('station',stations(ni),'lat',lat(ni),'lon',lon(ni),'e_rate',east(ni)/1000,'n_rate',north(ni)/1000, ...
                        'u_rate',up(ni)/1,'e_error',error_east(ni)/1,'n_error',error_north(ni)/1,'u_error',error_up(ni)/1 ...
                        ,'cov',zeros(size(stations,1)*3,size(stations,1)*3));
              end              
                
                
            case('default')

                [namelist]=textread(datafiles,'%s','delimiter','\n');

                % Check wether first line contains fields name

                %junk=list2str(namelist(1));
                %junkblk=find(junk==' ');

                if ~parameters
                    fprintf('No field info, assuming "lat lon vnorth veast vup errnorth erreast errup correlation station name" order \n');
                    tmp=[];

                    for ni=1:size(namelist,1)
                        junk=list2str(namelist(ni));
                        junkalph=isalpha(junk);
                        tmp=[tmp find(junkalph==1)];
                    end

                    tmp=unique(tmp);

                    for ni=1:size(namelist,1)
                        junk=list2str(namelist(ni));
                        GPS(ni).GPS_station={junk(tmp)};
                        junk(tmp)=[];junk=str2num(junk);junk=[junk repmat(NaN,1,9-size(junk,2))];
                        GPS(ni).date=[];
                        GPS(ni).dateYears=[];
                        GPS(ni).lat=junk(1);
                        GPS(ni).lon=junk(2);
                        GPS(ni).n_rate=junk(3)*unitsf;
                        GPS(ni).e_rate=junk(4)*unitsf;
                        GPS(ni).height_rate=junk(5)*unitsf;
                        GPS(ni).n_error=junk(6)*unitsf;
                        GPS(ni).e_error=junk(7)*unitsf;
                        GPS(ni).height_error=junk(8)*unitsf;
                        GPS(ni).corr=junk(9)*unitsf;
                        GPS(ni).Unit='m/yr';
                    end

                elseif parameters

                    % Get positions of fields
                    %stai=findstr('sta',junk);
                    lati=findstr('lat',parameters);loni=findstr('lon',parameters);vyi=findstr('vn',parameters);vxi=findstr('ve',parameters);vheighti=findstr('vu',parameters);
                    error_yi=findstr('sn',parameters);error_xi=findstr('se',parameters);error_heighti=findstr('su',parameters);corri=findstr('corr',parameters);

                    fieldvect=[lati;loni;vyi;vxi;vheighti;error_yi;error_xi;error_heighti;corri];
                    [fieldvects,isort]=sort(fieldvect);

                    % Get positions of station name

                    tmp=[];
                    for ni=1:size(namelist,1)
                        junk=list2str(namelist(ni));
                        junkalph=isalpha(junk);
                        tmp=[tmp find(junkalph==1)];
                    end
                    tmp=unique(tmp);
                    for ni=1:size(namelist,1)

                        junk=list2str(namelist(ni));
                        GPS(ni).GPS_station={junk(tmp)};
                        junk(tmp)=[];junk=str2num(junk);
                        GPS(ni).date=[];
                        GPS(ni).dateYears=[];
                        GPS(ni).lat=junk(find(isort==1));
                        GPS(ni).lon=junk(find(isort==2));
                        GPS(ni).n_rate=junk(find(isort==3))*unitsf;
                        GPS(ni).e_rate=junk(find(isort==4))*unitsf;
                        GPS(ni).height_rate=junk(find(isort==5))*unitsf;
                        GPS(ni).n_error=junk(find(isort==6))*unitsf;
                        GPS(ni).e_error=junk(find(isort==7))*unitsf;
                        GPS(ni).height_error=junk(find(isort==8))*unitsf;
                        GPS(ni).corr=junk(find(isort==9))*unitsf;
                        GPS(ni).Unit='m/yr';
                    end
                end
        end

    case('tmsr')

        if strcmp(source,'SOPAC')

            eval(sprintf('system(''ls %s > list1'')',datafiles));
            
        elseif strcmp(source,'Enrique')
            
            namelist = ListFiles(datafiles);
            
            for ni = 1:length(namelist)
                
                slashes             = max(strfind(namelist(ni).name,'/')) +1   ;  slashes(isempty(slashes)) = 1;
                GPS_Station         = namelist(ni).name(slashes:slashes+3)     ;
                fid                 = fopen(namelist(ni).name,'r')             ;
                C                   = textscan(fid,'%f%f%f%f','headerLines',1) ;
                GPS(ni).GPS_station = GPS_Station                              ;
                GPS(ni).dateYears  = C{1}                                     ;
                GPS(ni).e           = C{3}                                     ;
                GPS(ni).n           = C{2}                                     ;
                GPS(ni).u           = C{4}                                     ;
            end

        else

            str=pwd;
            inpath=str;

            count=1;
            %cd (path_GPS) % cd /RAID4/ngourmelen/GDs/
            %eval(sprintf('system(''ls %s > list'')',datafiles));
            %[namelist]=textread('list','%s','delimiter','\n');

            [namelist]=dir(sprintf('%s',datafiles));
            if ~namelist(1).isdir
                slashes=find(datafiles=='/');pwdir=datafiles(1:slashes(max(size(slashes))));
            else
                pwdir='./';
            end
            %namelist=getfields(namelist,'name','cell');
            namelist={namelist.name};namelist=namelist(:);
            num_GPS=size(namelist);
            num_GPS=num_GPS(1);

            for i=1:num_GPS

                name_GPS=namelist(i);
                name_GPS=[pwdir,char(name_GPS)];
                ss=find(name_GPS=='.');
                GPS_station=name_GPS(ss(1)-4:ss(1)-1);
                GPS_station={GPS_station}; % GPS station name

                text=textread(name_GPS,'%s','delimiter','\n');
                text=char(text);
                text(1,:)=[];

                num_mes=size(text,1);

                %%%%%%
                app=text(1,:);
                name_GPS;

                mes=str2num(app(11:73));

                lat=mes(1); lon=mes(3);
		
		if isempty(BoundingBox) 
   		  
		    for j=1:num_mes    
                        app=text(j,:);                        
			tt=find(app==' ');
                        dates(j,:)=app(1:tt(1)-1);
                        app_d1 = dates(j,:);
			app_d2 = ChangeDate(app_d1);                        
			app_d3 = datenum(app_d2,'yymmdd')/365.25;                        
			dates_gps(j,:)={app_d2};                        
			dates_years_gps(j,:)=app_d3;                        
			mes=str2num(app(11:73));

                        latitude(j,:)=mes(1);    
                        error_lat(j,:)=mes(2);                        
			longitude(j,:)=mes(3);
                        error_lon(j,:)=mes(4);                        
			height(j,:)=mes(5);
                        error_height(j,:)=mes(6);

                    end % for sulle misure
		    GPS_station
                    GPS(count).GPS_station=  GPS_station;                    
		    GPS(count).date=dates_gps;
                    GPS(count).dateYears=dates_years_gps;                    
		    GPS(count).lat=latitude*unitsf; 
                    GPS(count).lon=longitude*unitsf;
                    GPS(count).height=height*unitsf;
                    GPS(count).error_lat=error_lat*unitsf;                    
		    GPS(count).error_lon=error_lon*unitsf;                    
		    GPS(count).error_height=error_height*unitsf;
                    GPS(count).Unit='llm';
                    clear latitude error_lat longitude error_lon height error_height dates_gps dates_years_gps

                    if (plotgps)
                        figure
                        subplot(3,1,1);plot(GPS(count).dateYears,GPS(count).lat,'k*');title('Latitude') 
                        subplot(3,1,2);plot(GPS(count).dateYears,GPS(count).lon,'k*');title('Longitude')
                        subplot(3,1,3);plot(GPS(count).dateYears,GPS(count).height,'k*');title('height')
                    end
                    count=count+1;


                elseif (~exist('BoundingBox','var') || size(BoundingBox,2)==1 || (lon <= BoundingBox(2)) && (lon >= BoundingBox(1)) && (lat >= BoundingBox(3)) && (lat <= BoundingBox(4)))

                    for j=1:num_mes
                        app=text(j,:);
                        tt=find(app==' ');

                        dates(j,:)=app(1:tt(1)-1);

                        app_d1=dates(j,:);
                        app_d2=ChangeDate(app_d1);
                        app_d3=datenum(app_d2,'yymmdd')/365.25;

                        dates_gps(j,:)={app_d2};
                        dates_years_gps(j,:)=app_d3;
                        mes=str2num(app(11:73));

                        latitude(j,:)=mes(1);
                        error_lat(j,:)=mes(2);
                        longitude(j,:)=mes(3);
                        error_lon(j,:)=mes(4);
                        height(j,:)=mes(5);
                        error_height(j,:)=mes(6);

                    end % for sulle misure
                    GPS_station
                    GPS(count).GPS_station=  GPS_station;
                    GPS(count).date=dates_gps;
                    GPS(count).dateYears=dates_years_gps;
                    GPS(count).lat=latitude*unitsf;
                    GPS(count).lon=longitude*unitsf;
                    GPS(count).height=height*unitsf;
                    GPS(count).error_lat=error_lat*unitsf;
                    GPS(count).error_lon=error_lon*unitsf;
                    GPS(count).error_height=error_height*unitsf;
                    GPS(count).Unit='llm';
                    clear latitude error_lat longitude error_lon height error_height dates_gps dates_years_gps

                    if (plotgps)
                        figure
                        subplot(3,1,1);plot(GPS(count).dateYears,GPS(count).lat,'k*');title('Latitude')
                        subplot(3,1,2);plot(GPS(count).dateYears,GPS(count).lon,'k*');title('Longitude')
                        subplot(3,1,3);plot(GPS(count).dateYears,GPS(count).height,'k*');title('height')
                    end

                    count=count+1;

                end %if su lat e lon

            end %for sulle stazioni

        end

    case('matfile')
        GPS=load(datafiles);
        f=fieldnames(GPS) ;  fname={'GPS'};  eval([char(fname{1}) '= GPS.(f{1}) ;' ]) ;  
end

