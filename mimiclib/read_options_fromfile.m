function [optionsopt]=read_options_fromfile(fname,defaultopt);
%read_options_fromfile     - assigns file content to option structure
%
%usage: [optionsopt]=read_options_fromfile(fname,defaultopt);
%
%Input:   fname         text file containing options (*.min file)
%         defaultopt    option field names to be read
%         all fields 'off','on' are set to false,true
%
%Output:  optionsopt    new structure with options  
%
% Reads fields given in defaultopt from file. If defaultopt
% is not given reads the entire file.
% 
% TODO: Function needs to be completely revised. The first functionality may never be used and probably can be removed.
%       Now ignores everything following ### ROI_PAC OPTIONS ###. The should be replaced by function eliminating all roi_pac options using
%       and input list
%
%FA, Sep 2006    modified from process_defaultoptions.m
%FA, Aug 2007    made more stable. The input file is executed in matlab. All roi_pac options are ignored


optionsopt=[];
if ~exist(fname) return ;  end
if ~isempty(defaultopt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   f=fieldnames(defaultopt) ;  
   if ~exist('optionsopt','var')  optionsopt=[]  ; end
%   for i=1:length(f) optionsopt.(f{i}) = ReadMinFile(fname,f{i},'off') ; end   % FA 8/08  ReadMinFile is obsolete
%
%  now replace all 'off' by false
%
  f=fieldnames(optionsopt) ;  
   for i=1:length(f)
       if ~isstruct(optionsopt.(f{i})) &  strcmp('off',optionsopt.(f{i})) optionsopt.(f{i})=false ; end
   end
   % now replace all 'on' by true
   f=fieldnames(optionsopt) ;  
   for i=1:length(f)
       if ~isstruct(optionsopt.(f{i})) &  strcmp('on',optionsopt.(f{i})) optionsopt.(f{i})=true ; end
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
   % August 07 changes

   optionsopt=ReadKeywordfile(fname,'=');

   % code prior to Aug 07
   %fid=fopen(fname);  C = textscan(fid,'%s%s','delimiter','=','commentStyle','%') ; fclose(fid);
   %  names=deblank(C{1}');
   %  datas=deblank(C{2}');
   %  for i=1:length(names)
   %     %sprintf('%s=%s ;',names{i},data{i}) 
   %     %sprintf('optionsopt.%s=%s ;',names{i},data{i}) 
   %     eval(sprintf('%s=%s ;',names{i},datas{i})) ; 
   %     eval(sprintf('optionsopt.%s=%s ;',names{i},datas{i})) ; 
   %  end
end

