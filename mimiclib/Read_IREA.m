function [Fre] = Read_IREA(datafile,opt);

%
% Reads time series or rate map from IREA into a structure
%
%            [SBAS]=Read_IREA(datafile,opt)
%
%
% 'datafile'         : name of input file
%
%
% opt contains
%
% 'datesfile'    : name of the file containing the date list
%
% 'ResizeFactor' : Factor to resize the input matrice
%
%
%
%
% SBAS will be a structure, data will be in meters
%
%
% N. Gourmelen, March 2009
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultopt=struct(                                                ...
    'datesfile'                  ,        'off'      ,            ...
    'datatype'                   ,        'off'      ,            ...
    'ResizeFactor'               ,        'off'      )             ;


if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],[]); end ;
[opt]=process_defaultoptions(opt,defaultopt);  display(opt);
f=fieldnames(opt) ;
for i=1:length(f)
    eval([char(f{i}) '= opt.(f{i}) ;' ]) ;
end


%%%%%%%%%%%%%%%%%%%%%%%
% Check input variables

ipwd = find(datafile=='/') ;

if ipwd
    rootdir  = datafile(1:ipwd(end))                  ;
    datafile = datafile(ipwd(end)+1:length(datafile)) ;
else
    [status,rootdir] = unix('pwd') ;  rootdir = deblank(rootdir) ;
end

datafile = [rootdir,'/',datafile] ;

%%%%%%%%%%%%%%%%%%%%%%%%
% Setup some variables %
%%%%%%%%%%%%%%%%%%%%%%%%

if (~ResizeFactor)      ResizeFactor=1;                           end

[junk,junk,satbeam,sat]= extract_ProjectName(datafile) ; 
wavelength = extract_hardwired_satparameters(datafile, 'wavelength') ;

str=pwd;  inpath=str;

% Reads file with geographic infos

if exist([datafile,'.rsc'])
    rscinfo_ini    = ReadKeywordfile([datafile,'.rsc']) ;
    datatype       = rscinfo_ini.type                   ;
    rscinfo        = FilterKeyword(rscinfo_ini,'irea')         ;
    rscinfo.x_step = rscinfo.x_step/ResizeFactor ;
    rscinfo.y_step = rscinfo.y_step/ResizeFactor ;
else
    sprintf('WARNING:No parameter file named %s present',[datafile,'.rsc']) ;
end

if exist([datafile,'.msk'])
    maskfile  = [datafile,'.msk']   ;
    mask_fid  = fopen(maskfile,'r') ;
    [F,count] = fread(mask_fid,rscinfo.width*rscinfo.file_length,'float=>single') ;  fclose(mask_fid) ;
    mask      = reshape(F,[rscinfo.width,rscinfo.file_length]) ;
    mask      = rot90(mask) ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check whether rate or timeseries
if ~datatype datatype = 'generic' ; end

switch lower(datatype)
    case('interferograms')
        % Open file with interferograms
        if exist([datafile,'.dates'],'file')
            fid = fopen([datafile,'.dates']) ;  Call = textscan(fid,'%s','delimiter','\n') ;  fclose (fid) ;  Call = Call{1} ;
        else
            error('File %s does not exist',[datafile,'.dates'])
        end

        fid = fopen(datafile) ;  [F,count] = fread(fid,rscinfo.file_length*rscinfo.width*rscinfo_ini.NB_DATA,'float=>single') ;  fclose (fid) ;
        Fretmp     = reshape(F,[rscinfo.width,rscinfo.file_length,rscinfo_ini.NB_DATA]) ;  clear F ;

        for ni=2:length(Call)
            Call_tmp              = Call{ni}                  ;
            Fre(ni-1).data(:,:)   = rot90(Fretmp(:,:,(ni-1))) ;
            Fre(ni-1).width       = size(Fretmp,2)            ;
            Fre(ni-1).file_length = size(Fretmp,1)            ;
            Fre(ni-1).date1       = [Call_tmp(5:8),Call_tmp(3:4),Call_tmp(1:2)]       ;
            Fre(ni-1).date2       = [Call_tmp(17:20),Call_tmp(15:16),Call_tmp(13:14)] ;
            Fre(ni-1).t1          = date2j(str2num(Fre(ni-1).date1(3:4)),str2num(Fre(ni-1).date1(5:6)), ...
                str2num(Fre(ni-1).date1(7:8)))          ;
            Fre(ni-1).t2          = date2j(str2num(Fre(ni-1).date2(3:4)),str2num(Fre(ni-1).date2(5:6)), ...
                str2num(Fre(ni-1).date2(7:8)))          ;
            Fre(ni-1).delt        = Fre(ni-1).t2 - Fre(ni-1).t1   ;
            Fre(ni-1).sat         = satbeam 						  ;
			Fre(ni-1).sat_height  = 793882.7391                   ;
            Fre(ni-1).Unit        = rscinfo_ini.Unit              ;
            tmpFre(ni-1)          = add_struct(Fre(ni-1),rscinfo) ;
        end
        Fre = tmpFre ;
    case('tmsr')
        % Open date file
        if exist([datafile,'.dates'],'file')
            fid = fopen([datafile,'.dates']) ;  dates = fread(fid,'double') ;  fclose (fid) ;
        else
            error('File %s does not exist',[datafile,'.dates'])
        end

        fid = fopen(datafile) ;  [F,count] = fread(fid,rscinfo.file_length*rscinfo.width*rscinfo_ini.NB_DATA,'float=>single') ;  fclose (fid) ;
        Fretmp     = reshape(F,[rscinfo.width,rscinfo.file_length,rscinfo_ini.NB_DATA]) ;  clear F ;
        fact       = convert_unit('m',rscinfo_ini.Unit) ;
        for ni=1:rscinfo_ini.NB_DATA
            Fre(ni).data(:,:)   = rot90(Fretmp(:,:,(ni))) * fact ;
            Fre(ni).width       = size(Fretmp,2)                 ;
            Fre(ni).file_length = size(Fretmp,1)                 ;
            Fre(ni).date        = y2yymmdd(dates(ni))            ;   % y2yymmdd maybe not that precise!
            Fre(ni).dateYears   = dates(ni)                      ;
            if ni==1 to = date2j(str2num(Fre(ni).date(1:1)),str2num(Fre(ni).date(3:4)), ...
                str2num(Fre(ni).date(5:6)))   ;  end
            Fre(ni).t          = date2j(str2num(Fre(ni).date(1:1)),str2num(Fre(ni).date(3:4)), ...
                str2num(Fre(ni).date(5:6)))    - to      ;
			Fre(ni).sat         = satbeam                           ;
            Fre(ni).sat_height  = 793882.7391         ;
            tmpFre(ni)          = add_struct(Fre(ni),rscinfo) ;
            tmpFre(ni).Unit     = 'm'                         ;
        end
        Fre = tmpFre ;
    case('rate')
        fid         = fopen(datafile,'r') ;             [F,count]       = fread(fid,'float=>single') ;  fclose(fid);
        fact        = convert_unit('m/yr',rscinfo_ini.Unit) ;
        Fre.data    = reshape(F,[rscinfo.width,rscinfo.file_length]) * fact;
        Fre.data    = rot90(Fre.data)     ;
        Fre.width   = size(Fre.data,2)    ;             Fre.file_length = size(Fre.data,1) ;
        Fre.Unit    = 'm/yr'              ;
		Fre.sat         = satbeam                           ;
        Fre         = add_struct(Fre,rscinfo) ;
        Fre.time_range = [yymmdd2y(Fre.date1) yymmdd2y(Fre.date2)] ;
    case('generic')
        fid         = fopen(datafile,'r') ;             [F,count]       = fread(fid,'float=>single') ;  fclose(fid);
        Fre.data    = reshape(F,[rscinfo.width,rscinfo.file_length]) ;
        Fre.data    = rot90(Fre.data)     ;
        Fre.width   = size(Fre.data,2)    ;             Fre.file_length = size(Fre.data,1) ;
        Fre.Unit    = rscinfo.Unit        ;
		Fre.sat         = satbeam                           ;
        Fre         = add_struct(Fre,rscinfo) ;
end




    
    
