function [Quakes,Seismi] = Read_Seismicity(datafile,opt);

%
% reads seismicity file
%
% [Quakes,Seismi] = Read_Seismicity(struct('filetype','rates','datafiles','/RAID1/ngourmelen/GPSdata/YuccaProfile.dat','units','mm'));
%
%      Output is Quakes=[lon,lat,depth,mag] or Seismi=structure with fields name
%
%
%   datafile     :   File containing Seismicity data
%
%
% opt is:
%
%
% parameters    :   Vector with keywords: [lat lon dep mag yrs yr month day]
%
% source	    :   'local', 'anss_global', 'anss', 'harvard', 'engdahl', 'LinReloc', ...
%
% location   :   [lat1 lat2 lat3 lat4 ...;lon1 lon2 lon3 lon4 ...] are the geographic coordinates of the study area polygon
%
% startT and endT    :   start and end Time of the form 'yyyymmdd'
%
% startM and endM    :   start and end Magnitude
%
% startD and endD    :   start and end Depth
%
%
% N. Gourmelen, September 2006
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultopt=struct(                                               ...
    'parameters'               ,        'off'           ,        ...
    'source'                   ,        'off'	    	,	     ...
    'startT'                   ,        'off'           ,        ...
    'endT'                     ,        'off'           ,        ...
    'startM'                   ,        'off'           ,        ...
    'endM'                     ,        'off'           ,        ...
    'startD'                   ,        'off'           ,        ...
    'endD'                     ,        'off'           ,        ...
    'location'                 ,        'off'           ,        ...
    'plotseism'                ,        'off'      )             ;

%
%  read *min file and process  options (selectpairsopt)
%

if ~isstruct(opt)
    [S]=ReadKeywordfile(opt,'=') ;[opt]=process_defaultoptions(S.readseismiopt,defaultopt);
else
    [opt]=process_defaultoptions(opt,defaultopt);
end

disp (opt)

f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end



%%%%%%%%%%%%%%%%%%%
% Check variables %
%%%%%%%%%%%%%%%%%%%

switch lower(source)
    
    case('local')
        
        if (datafile & parameters)    
            
        else            exit ('Error, no datafile specified');  end
    
    case('anss_global')
        
        matfile  ='/RAID1/ngourmelen/SEISMIdata/ANSSglobal.mat';
        filtfile ='/RAID1/ngourmelen/SEISMIdata/ANSSglobal.dat.filt';
        datfile  ='/RAID1/ngourmelen/SEISMIdata/ANSSglobal.dat';
        
        %if (exist(matfile)==2)
            
        %    sprintf('Load %s ',matfile); 
        %    load (matfile);
            
        if (exist(filtfile)==2)
            
            fid=fopen(filtfile,'r');C=textscan(fid,'%s %s %f %f %f %f %f');fclose(fid);
            
            Seismi.years=C{7};Seismi.mag=C{6};Seismi.depth=C{5};Seismi.lat=C{3};Seismi.lon=C{4};Seismi.date=C{1};Seismi.time=C{2};
            
            %save /RAID1/ngourmelen/SEISMIdata/ANSSglobal Seismi
        
        else
            
            exit('Error, no data file at /RAID1/ngourmelen/SEISMIdata/');
            
        end
        
    case{'anss_readable'}
        fid=fopen(datafile); C=textscan(fid,'%s %s %f %f %f %f %s %s %s %s %s %s %s'); 
        
        Seismi.years=C{7};Seismi.mag=[C{6}'];Seismi.depth=[C{5}'];Seismi.lat=[C{3}'];Seismi.lon=[C{4}'];Seismi.date=[C{1}'];Seismi.time=[C{2}'];
        
        Quakes=[long',lat',depth',mag'];
        % Notes:
        %We need to read from the file according to the location (character) which is unchanged over the entire file.
        %Also, how to specify header lines with textscan ?
        %
        %One way to to this would be:
        %C=textscan(fid,'%s','delimiter','\n')
        %tmp1= char(C{1}');
        %lat   = str2num(tmp1(:,25:30));
        %
        %but is this efficient ?
        %
        %Some information may be given my matlab's importwizard.
        
    case{'anss'}
        fid=fopen(datafile); C=textscan(fid,'%f%f%f%s%f%f%s%s%s%s%s%s%s%s%s')
        
        Seismi.years=C{7};Seismi.mag=[C{6}']*0+1;Seismi.depth=str2num(tmp2(:,9:13));Seismi.lat=str2num(tmp1(:,9:15));
        Seismi.lon=[C{5}'];Seismi.date=[C{2}'];Seismi.time=[C{2}'];
        year  = [C{1}'];
        month = [C{2}'];
        day   = [C{3}'];
        tmp1= char(C{4}');
        tmp2= char(C{5}');         % This does not work. Don't know why. Not necessary to fix if we use anss_readable
        mag   = [C{6}'];
        lat   = str2num(tmp1(:,9:15));
        %long  = str2num(tmp2(:,1:7));
        depth = str2num(tmp2(:,9:13));
        Quakes=[long,lat,depth,mag];

            
    case{'INGV-CT-FocMecs'}
        fid=fopen(datafile); C=textscan(fid,'%d %f %s %f %f %f %s %s %s %s %s %s','headerLines',2); fclose(fid);
        
        Seismi.years=C{7};Seismi.mag=[C{6}']*0+1;Seismi.depth=[C{6}'];Seismi.lat=[C{4}'];Seismi.lon=[C{5}'];Seismi.date=[C{2}'];Seismi.time=[C{3}'];
      

    case('sod_iris')
        
        matfile='/RAID1/ngourmelen/SEISMIdata/SOD_IRIS.mat';
        datfile='/RAID1/ngourmelen/SEISMIdata/SOD_IRIS.dat';
        
        if (exist(matfile)==2)
            
            load matfile;
            
        elseif (exist(datfile)==2)

            indata=textread(datfile,'%s','delimiter','\n');
            
            for ni=1:length(indata)
                
                
            end
                
            save /RAID1/ngourmelen/SEISMIdata/ANSSglobal.mat Seismi;
        
        else
            
            exit('Error, no data file at /RAID1/ngourmelen/SEISMIdata/');
            
        end
     

            
    case('engdahl')
        
        indata=textread(datafile,'%s','delimiter','\n');
        
        for ni=1:length(indata)
            
            ff=indata{ni};
            Seismi.date(ni)  =datenum(str2num(ff(12:15)),str2num(ff(17:18)),str2num(ff(20:21)))/365.25;
            Seismi.mag(ni)   =str2num(ff(67:69));
            Seismi.depth(ni) =str2num(ff(53:57));
            Seismi.lat(ni)   =str2num(ff(37:43));
            Seismi.lon(ni)   =str2num(ff(44:51));
            
        end
    
    case('neic')
        
        C = textscan(fid,'%f%f%f%f%u64%f','Delimiter',',','EmptyValue',-Inf);
        [indata]=textread(datafile,'%s','delimiter','\n');
        
    case('linreloc')
        
        fid = fopen(datafile) ;
        C   = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s','Delimiter',',','EmptyValue',-Inf) ;
        Seismi.years =  C{1} ;
        Seismi.month =  C{2} ;
        Seismi.day   =  C{3} ;
        Seismi.lat   =  C{8} ;
        Seismi.lon   =  C{9} ;
        Seismi.depth = C{10} ;
        Seismi.mag   = C{11} ;
        
end

in=[];

if startT    stTin=find(Seismi.years<datenum(str2num(startT(1:4)),str2num(startT(5:6)),str2num(startT(7:8)))/365.25);       in=unique([in;stTin]);   end
    
if endT      enTin=find(Seismi.years>datenum(str2num(endT(1:4)),str2num(endT(5:6)),str2num(endT(7:8)))/365.25);             in=unique([in;enTin]);   end

if location  locin=find(inpolygon(Seismi.lat,Seismi.lon,location(1,:),location(2,:))==1);                                   in=unique([in;locin]);   end

if startM    stMin=find(Seismi.mag    <  startM);                                                                           in=unique([in;stMin]);   end

if endM      enMin=find(Seismi.mag    >    endM);                                                                           in=unique([in;enMin]);   end

if startD    stDin=find(Seismi.depth  <  startD);                                                                           in=unique([in;stDin]);   end

if endD      enDin=find(Seismi.depth  >    endD);                                                                           in=unique([in;enDin]);   end


Fields = fieldnames(Seismi) ;

for ni=1:length(Fields)
    
    tmp_data = Seismi.(genvarname(Fields{ni})) ;
    Seismi.(genvarname(Fields{ni})) = tmp_data(in) ;
    
end

%Seismi.years=Seismi.years(in);Seismi.lon=Seismi.lon(in);Seismi.lat=Seismi.lat(in);Seismi.date=Seismi.date(in);Seismi.time=Seismi.time(in);Seismi.mag=Seismi.mag(in);Seismi.depth=Seismi.depth(in);


Quakes=[Seismi.lon,Seismi.lat,Seismi.depth,Seismi.mag];
    






