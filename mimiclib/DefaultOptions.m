function [S]=DefaultOptions(fname,defaultopt);
%DefaultOptions     - assigns dir names according to RSMAS naming convention
%
%usage:   [S] = DefaultOptions(fname,options);
%
%Input:   fname         text file containing options (*.min file)
%
%Output:  S    new structure with options  
%
% 
% TODO: needs revamping so that it does all the option processing for tssar, SelectPairs and geodmod
%       read_options_from_file should be eliminated and al;so the command option processing should be done here
%
%FA, Oct 2009    

 S                     = [];
 [path,name,ext]       = fileparts(fname);
 if isempty(ext) fname = fullfile(path,[name '.min']); end    % FA 10/2009: should check whether *template or *min file exist and use fname

 S                     = ReadKeywordfile(fname,'=');

 if ~exist('ProcessDir','var')       ProcessDir = getenv('PROCESSDIR'); end 
 if ~exist('TssarDir','var')         TssarDir   = getenv('TSSARDIR')  ; end 
 if  isfield(S,'dir_out_parent')     ProcessDir = S.dir_out_parent    ; end
 if  isfield(S,'dir_out_parent')     TssarDir   = S.dir_out_parent    ; end

 if ~isfield(S.inputdataopt,'interf') S.inputdataopt.interf = [ ProcessDir filesep 'process' filesep name filesep 'DONE/P*/geo_*-*_*_*.unw' ]      ; end
 if ~isfield(S.inputdataopt,'bllist') S.inputdataopt.bllist = [ TssarDir   filesep name filesep 'DONE/bl_list.txt' ]             ; end
 if ~isfield(S.inputdataopt,'topo')   if isfield(S,'DEMg') S.inputdataopt.topo = S.DEMg; else S.inputdataopt.topo = S.DEM; end ; end

