function []=SaveProject(fname,location)
%SaveProject         - copies data, dem files, etc into dedicated directory and adjusts *min file 
%
%  usage:  SaveProject(fname,location,flag_relative)
%
%     e.g. SaveProject('Kil2005-2007.min')
%          SaveProject('Kil2005-2007.min','pd')
%          SaveProject('Kil2005-2007.min','ex')
%          SaveProject('Kil2005-2007.min','rel')
%
%          the new *min file is called Kil2005-2007_cp.min  or Kil2005-2007_ex.min
%
%   Options:
%
%          fname:         geodmod control file  (*.min)
%          location:      location where data are copied  ('rel','pd','ex')  [default 'pd' ]
%                         'rel'   project_data (relative path)
%                         'pd'    environment variable GEODMOD_PROJECTDATA (absolute path)
%                         'ex'    INSARLABHOME/testdata/geodmod (only for exceptional cases)
%
%          TODO: Need an option to copy only data but not DEM, etc                                                     
%          TODO: for quake files it should check for the size and change *min file only if size is below a certain limit
%                obviously, there need to be 'force' option
%
% March 2007 Falk Amelung
%

if nargin==0  help SaveProject; end
[opt]=read_options_fromfile(fname,[]); f=fieldnames(opt); for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end ; 

if ~exist('location','var')       location='project_data' ; end
flag_corfile_done = false ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                  %
% First set up destination directory and filenames %
%                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pathstr, fname_stem, ext] = fileparts(fname);

[s,GEODMOD_PROJECTDATA]  = unix('echo $GEODMOD_PROJECTDATA'); GEODMOD_PROJECTDATA = [deblank(GEODMOD_PROJECTDATA) filesep];
[s,INSARLABHOME       ]  = unix('echo $INSARLABHOME       '); INSARLABHOME        = [deblank(INSARLABHOME       ) filesep];

if strcmp(location,'rel') || strcmp(location,'relative')
      suffix = '_rel' ;
      rootdir= 'project_data' ;
elseif strcmp(location,'pd') || strcmp(location,'project_data')
      suffix = '_cp' ;
      rootdir= GEODMOD_PROJECTDATA ;
elseif strcmp(location,'ex') || strcmp(location,'example')
      suffix = '_ex' ;
      rootdir=  fullfile( INSARLABHOME,'testdata','geodmod' ) ;
else
      error('not adopted for location %s',location) 
end

newpath= fullfile (rootdir,fname_stem) ;
if ~exist(newpath,'dir')  mkdir(newpath) ; end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  list of strings for which line in *.min file needs to be changed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savelist_SAR      = { 'makesaropt.dir_inlist{1}'                   ...
                      'makesaropt.dir_inlist{2}'                   ...
                      'makesaropt.dir_inlist{3}'                   ...
                      'makesaropt.dir_inlist{4}'                   ...
                      'makesaropt.dir_inlist{5}'                   ...
                      'makesaropt.dir_inlist{6}'                   ...
                      'makesaropt.dir_inlist{7}'                   ...
                    };
savelist_GPS      = { 'prepareopt.gpsfile'          ...
                      'prepareopt.GPSfile'          ...
                      'makedatasetopt.gpsfile'           ...
                      'makedatasetopt.GPSfile'           ...
                    };
savelist_dem      = { 'prepareopt.demfile'          ...
                      'prepareopt.ShadedRelieffile' ...
                    };
savelist_lines    = { 'prepareopt.linefile'         ...
                    };
savelist_quakes   = { 'prepareopt.quakefile'        ...
                    };
savelist_all      = [savelist_SAR savelist_GPS savelist_dem savelist_lines savelist_quakes];

specialtreat_list = {  'makesaropt.roipac2igramopt.pattern'        ...
                    };

%
% open output file file (with suffix _cp or _ex )
%

[pathstr, name, ext] = fileparts(fname);
nname = fullfile(pathstr,[name suffix ext]);
fid1=fopen(nname,'w');

%
% read *.min file into 2-column cell array and loop over all lines
%

fid2=fopen(fname);  C = textscan(fid2,'%s%s','delimiter','=') ; fclose(fid2); clear fid2 ;

for i=1:length(C{1})
    if isempty(cell2mat(strfind(savelist_all,deblank(char(C{1}(i))))))                 % empty if no match with contents of savelist
        value =  char( C{2}(i)) ;                                                  % copy content of orginal file
    else                                                                           % rename and copy fields found in savelist
        logmessage(sprintf('Working on: %s',char( C{2}(i))))
        oldfname  =char( C{2}(i));
        oldfname  =eval(char( C{1}(i)));
        [pathstr,name,ext] = fileparts(oldfname);
        if     ~isempty(cell2mat(strfind(savelist_SAR  ,deblank(char(C{1}(i))))))     % empty if no match with contents of savelist_SAR  
                          ndir=fullfile(newpath,'Data') ;  if ~exist(ndir,'dir')  mkdir(ndir) ; end 


                %
                % copy geo*unw files
                %
                          logmessage(['copying data to ' ndir])
                          filematch=fullfile(eval(deblank(char(C{1}(i)))),makesaropt.roipac2igramopt.pattern);
                          [junk,junk,suffix]=extract_name_from_SOdir(filematch);
                          ndir=fullfile(newpath,'Data',suffix) ;  if ~exist(ndir,'dir')  mkdir(ndir) ; end 
                          newfname=ndir;
                          commandstring = sprintf('cp %s %s',filematch,ndir);      [status,value]=unix(commandstring);
                          rsc_commandstring=strrep(commandstring,'unw','unw.rsc'); [status,value]=unix(rsc_commandstring);
                %
                % copy geo*cor files
                %
                          logmessage('copying coherence files')
                          commandstring=strrep(commandstring,'unw','cor');         [status,value]=unix(commandstring);
                          rsc_commandstring=strrep(commandstring,'cor','cor.rsc'); [status,value]=unix(rsc_commandstring);

                %
                % copy mean*cor file or similar
                %
                          if isfield(makesaropt.roipac2igramopt,'corfile') &&  ~ flag_corfile_done
                             logmessage('copying special (mean) coherence file')
                             filematch=fullfile(eval(deblank(char(C{1}(i)))),roipac2igramopt.corfile);
                             commandstring = sprintf('cp %s %s',filematch,ndir);      [status,value]=unix(commandstring);
                             rsc_commandstring=strrep(commandstring,'cor','cor.rsc'); [status,value]=unix(rsc_commandstring);
                             flag_corfile_done=true;
                          end
                     
        elseif ~isempty(cell2mat(strfind(savelist_GPS   ,deblank(char(C{1}(i))))))     
                ndir=fullfile(newpath,'Data') ;  if ~exist(ndir,'dir')  mkdir(ndir) ; end 
                newfname= fullfile(ndir,[name  ext]) ;
                logmessage(['copying GPSdata to ' newfname])
                copyfile( oldfname,         newfname        );
        elseif ~isempty(cell2mat(strfind(savelist_dem   ,deblank(char(C{1}(i))))))     
                ndir=fullfile(newpath,'Dem') ;  if ~exist(ndir,'dir')  mkdir(ndir) ; end 
                newfname= fullfile(ndir,[name  ext]) ;
                logmessage(['copying demfile to ' newfname])
                copyfile( oldfname,         newfname        );
                copyfile([oldfname '.rsc'],[newfname '.rsc']);
        elseif ~isempty(cell2mat(strfind(savelist_lines ,deblank(char(C{1}(i))))))     
                ndir=fullfile(newpath,'Lines') ;  if ~exist(ndir,'dir')  mkdir(ndir) ; end 
                newfname= fullfile(ndir,[name  ext]) ;
                logmessage(['copying linefile to ' newfname])
                copyfile(oldfname,newfname);
        elseif ~isempty(cell2mat(strfind(savelist_quakes,deblank(char(C{1}(i))))))     
                ndir=fullfile(newpath,'Quakes') ;  if ~exist(ndir,'dir')  mkdir(ndir) ; end 
                newfname= fullfile(ndir,[name  ext]) ;
                logmessage(['copying quakefile to ' newfname])
                copyfile(oldfname,newfname);
        end

        value =  [ '''' newfname '''' ];
    end
    if strfind(deblank(char(C{1}(i))),'makesaropt.roipac2igramopt.pattern')           % empty if no match with contents of savelist
        [pathstr,name,ext] = fileparts(makesaropt.roipac2igramopt.pattern);
        newname=fullfile('',[name  ext]) ;
        value =  [ '''' newname '''' ];
    end
       
    if ~isempty(char(C{2}(i))) 
        str=[ char(C{1}(i)) '=  ' value ] ;
    else 
       str=[ char(C{1}(i))  ] ;             % e.g. comment lines
    end
    fprintf(fid1,'%s\n',str);
end 
