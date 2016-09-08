function [InSARout,opt]=CheckDownloadData(prepareopt,makesaropt);
%CheckDownloadData     -  Checks whether data dirs exist. Download if necessary
%
% Downloads *tgz file from url if Dem and Data are not found.
%
% Convention: if Dem is in  /RAID6/insar_lab/testdata_geodmod/Wells/Dem/wells_srtm3.dem
%             and Data in   /RAID6/insar_lab/testdata_geodmod/Wells/Data/EnvD2
%             then  Wells.tgz is downloaded and untared in /RAID6/insar_lab/testdata_geodmod/
%              
%
% remote server is specified in prepareopt.url
%
% FA, Oct 2008
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultopt=struct(                          ...
    'url'                ,  'http://www.rsmas.miami.edu/personal/famelung/geodmod'   ,   ...
    'plotres'            ,  'off'   );

[prepareopt]=process_defaultoptions(prepareopt,defaultopt);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if dem file and data dir exist. Return if they do.                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if exist(prepareopt.demfile,'file') && (~isfield(makesaropt,'dir_inlist') || ( ischar(makesaropt.dir_inlist) && exist(makesaropt.dir_inlist,'dir')  ||  iscell(makesaropt.dir_inlist) && exist(makesaropt.dir_inlist{end},'dir')))
        logmessage('Dem file and Data directories found'); return
    end
%
        dem_dir         = prepareopt.demfile(1:findstr('Dem', prepareopt.demfile)+2);
        dem_file        = prepareopt.demfile(  findstr('Dem', prepareopt.demfile)+4 : end);
        dem_dir_parent  = fileparts(dem_dir); 

    if isfield(makesaropt,'dir_inlist')                                                                   %dir_inlist may not be given (e.g. Forward_mogi)
        data_dir        = makesaropt.dir_inlist{1}(1:findstr('Data', makesaropt.dir_inlist{1})+3);
        data_file       = makesaropt.dir_inlist{1}(  findstr('Data', makesaropt.dir_inlist{1})+5 : end);
        data_dir_parent = fileparts(data_dir); 
    else
        data_dir_parent = dem_dir_parent; 
        data_dir        = dem_dir;
    end

       if ~strcmp(data_dir_parent, dem_dir_parent)  errordlg('Different parent directories for Dem and Data'); error('user error -- exiting'); end
        
        [projdir_parent,proj_name] = fileparts(dem_dir_parent);
        
        remote_datafile            = [prepareopt.url filesep proj_name '.tgz'];
        logmessage('##########################################################')
        logmessage(['Downloading ' remote_datafile ' ....']); 
        logmessage('##########################################################')
        
        untar( remote_datafile, projdir_parent);
    end
