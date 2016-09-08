function [igram,datelist,N_timevector,t]=LoadData(datafile,fname)


%  LoadData  -  Load data into igram structure ...
%
%  usage:  [igram,datelist,N_timevector]=LoadData(LoadOpt);
%
%          'datafile'    InFiles
%
%   LoadOpt:
%
%          'maskfile'      Used for IREA loading only
%          'datesfile'     Used for IREA loading only
%          'origin'        Type of file: 'roi_pac', 'irea', 'doris', ...
%          'purpose'       'geodmod, 'SBAS', ...
%          'subset'        Subset the Interferograms using resize_igram.m
%          'loadAmplitude' Loads the amplitude along with phase for roi_pac case
%          'dataMean'      returns only the mean of the data (useful for reading geo_incidence.unw)
%          'ReadInteg'     Use Integer format to store the interferograms
%          'readformat'    'single', ... Use Single, ... format to store the interferograms
%          'DryRun'        option to only check a few parameters (so far only "Unit")
%
%
%  Part of MILK (Miami InSAR Lab Kernel) 
%
%  N Gourmelen, November 2005
%   revision:
%           N Gourmelen, August 2007
%


%% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultopt=struct(                                    ...
    'origin'         ,        'off'      ,            ...
    'purpose'        ,        'off'      ,            ...
    'subset'         ,        'off'      ,            ...
    'loadAmplitude'  ,        'off'      ,            ...
    'dataMean'       ,        'off'      ,            ...
    'DryRun'         ,        'off'      ,            ...
    'ReadInteg'      ,        'off'      ,            ...
    'readFormat'     ,        'double'   )            ;


if nargin==0 fname=[mfilename '.min']; end

if ~isstruct(fname)
    [S]   = ReadKeywordfile(fname,'=') ;  [opt] = process_defaultoptions(S.readseismiopt,defaultopt) ;
else
    [opt] = process_defaultoptions(fname,defaultopt) ;
end

f = fieldnames(opt) ; for ni = 1:length(f) eval([char(f{ni}) '= opt.(f{ni}) ;' ]) ;  end ;  

%% Get File name
infiles = ListFiles(datafile) ;

%% Variables definition 
if strcmp(datafile(end-2:end),'mat')  origin = 'matfile' ;  end
test='junk' ;  if (~origin);  origin = 'roi_pac' ;  end ;  if (~purpose);  purpose = 'geodmod' ;  end

%% Load The Interfero %
switch lower(origin)
    case('matfile')
        Stmp = load(datafile) ;  igram = Stmp.(genvarname(list2str(fieldnames(Stmp))));
        if isstruct(subset)  
            for ni=1:length(igram)
                opt.subset = subset ;  igram(ni) = resize_igram(igram(ni),opt) ;  
            end
        end
        if DryRun  igram = igram(1).Unit;  return ;  end ;
        %if isstruct(subset)  opt.subset = subset ;  igram(ni) = resize_igram(igram(ni),opt) ;  end
        if isfield(igram(1),'Unit')  opt.Unit = igram(1).Unit ;  end
        
    case('irea')
        opt = [];  opt = struct('datafile',datafile);
        [igram] = Read_IREA(datafile,opt);
        if DryRun  igram = igram(1).Unit;  return ;  end ;
        if isstruct(subset)  opt.subset = subset ;  igram(ni) = resize_igram(igram(ni),opt) ;  end

    case('doris')
        
    case('roi_pac')
        % Get minimum lenght YMAX0 of the interfero %
        YMIN0 = 1e32 ;
        for ni=1:length(infiles)
            RscInfo = ReadKeywordfile([infiles(ni).name,'.rsc']) ;
            if RscInfo.FILE_LENGTH < YMIN0        YMIN0 = RscInfo.FILE_LENGTH ;           else          end
        end
        
        if DryRun ;  infiles = infiles(1) ;  end ;
        
        for ni=1:length(infiles)
            file     = infiles(ni).name                                                                           ;
            str      = sprintf('Loading file %d %s',ni, file) ; logmessage(str)
            [Infos]  = ReadKeywordfile([file,'.rsc']);  [Infos] = FilterKeyword(Infos,origin,purpose)             ;

            [a,p]    = readfile(file,struct('precision',readFormat))       ; p(find(a==0))        = NaN ; a = a(1:YMIN0,:) ;  p = p(1:YMIN0,:)   ;
            igramtmp = Infos ; igramtmp.data = p ; 
            
            if loadAmplitude ;  
                igramtmp.amp = a ;
                clear a;
            end ;
            if dataMean igramtmp.dataMean=[]; end                                                                 % FA 1/2011: need to be initialized so that resize_igram assignement works
            %if strcmp(readFormat,'single')  
            %    igramtmp.data        =   single(p) ;
            %else
            %    igramtmp.data        =   p ;
            %end
           
            junk     = fieldnames(igramtmp) ; junkl=[1:length(junk)]                                              ; %FA 1/2011: the 3 lines seem unncessary expect Unit='Radian'
            junkl(1) =  length(junk)        ; junkl(length(junk))  =   1 ; igramtmp = orderfields(igramtmp,junkl) ;
            test = file((length(file)-2):length(file));
            
            igramtmp.Unit = 'radian'; 

            if strcmp(purpose,'sbas')

                slash=find(file==filesep);baselinerscfile=strcat(file(1:slash(size(slash,2))),igramtmp.date1,'_',igramtmp.date2,'_baseline.rsc');

                if exist(baselinerscfile)
                    [Infos]  = ReadKeywordfile(baselinerscfile)                ;
                    [Infos]  = FilterKeyword(Infos,'roi_pac_baseline',purpose) ;
                    Infos.file_length = YMIN0 ;
                    igramtmp = add_struct(igramtmp,Infos)                      ;
                else
                    error('No baseline file %s',baselinerscfile)
                end
            end
            
            [junk,junk,satbeam,sat] = extract_ProjectName(infiles(ni).name) ;  igramtmp.sat = satbeam ;  
            
			if isstruct(subset)  
                opt.subset=subset;  
                opt.field ={'data'}
                if loadAmplitude
                    opt.field = {'data' 'amp'};
                end
                igram(ni) = resize_igram(igramtmp,opt) ;
               % if loadAmplitude
               %    tmpopt       = opt;
               %    tmpopt.field = 'amp';
               %    igram(ni)    = resize_igram(igramtmp,tmpopt) ;
               % end                   
            else
                igram(ni) = igramtmp ;  
            end
            
            if dataMean; 
                igram(ni).dataMean=nanmean(nanmean(igram(ni).data));   clear igram.data; 
            end
            
            if (ni==1);  igram(length(infiles)) = igram(1) ;  end  % Use to initialize structure
            
        end        
        
        if DryRun ;  igram = igram(1).Unit ;  datelist = 'junk' ;  N_timevector = 'junk' ; t = 'junk' ;  end ;
        
    case('gamma')
        
        
        
    case('dem')
        
        [junk,dem_data,Infos]=readfile(datafile)    ;  %[Infos]=ReadKeywordfile([datafile,'.rsc']);  
        [Infos]=FilterKeyword(Infos,origin,purpose) ;  igram = Infos ;  igram.data=dem_data ;
        if isstruct(subset)  opt.subset=subset      ;  igram = resize_igram(igram,opt)      ;  end
end

if DryRun ;  return ;  end ;

%%  convert phase into unit desired
if ~isfield(opt,'Unit')    opt.Unit = 'radian'        ;  end
if (strcmp(test,'cor'))    opt.Unit =  'unity'        ;  end
if (strcmp(origin,'irea')) opt.Unit =  igram(1).Unit  ;  end
if (strcmp(origin,'dem'))  opt.Unit =  igram(1).Unit  ;  end
    
fac = convert_unit(opt.Unit,igram(1).Unit) ;
for ni=1:length(igram)  igram(ni).data = igram(ni).data*fac ;  igram(ni).Unit = opt.Unit ;  end

%% sort interferograms 
if isfield(igram(1),'t1')
    [igram,N_igrams] = SortIgram(igram) ;
    
    % set time to first image to zero and create datelist %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ni=1:N_igrams  t(ni,:) = [igram(ni).t1 igram(ni).t2] ;  end
    firsttime = min(min(t)) ;
    for ni=1:N_igrams
        igram(ni).t1 = igram(ni).t1 - firsttime    ;  igram(ni).t2 = igram(ni).t2 - firsttime ;   
        t(ni,:)      = [igram(ni).t1 igram(ni).t2] ;  
    end
    datelist = unique(t) ;  N_timevector = length(datelist) ;
end
if isfield(igram(1),'t')
    [igram,N_igrams] = SortIgram(igram) ;
    
    % set time to first image to zero and create datelist %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ni=1:N_igrams  t(ni,:) = [igram(ni).t] ;  end
    firsttime = min(t) ;
    for ni=1:N_igrams
        igram(ni).t = igram(ni).t - firsttime ;  t(ni,:)      = [igram(ni).t] ;  
    end
    datelist = unique(t) ;  N_timevector = length(datelist) ;
end
