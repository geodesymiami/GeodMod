function []=gclean(what, arg2);
%gclean     -  removes files from last geodmod run
%
%usage:  gclean igram
%        gclean igram only       
%
%   For the following options the following files are removed:
%
%   gclean all           removes all
%   gclean model         model
%   gclean igram         igram_*
%   gclean motion        motion_*
%   gclean timeseries    timeseries_*
%   gclean Qt            motion*Qt*
%   gclean Grid          motion*Grid*
%   gclean Profile       motion*Profile*
%   gclean threed        motion_threedfield
%   gclean GPS           *GPS*
%   gclean basemap       basemap
%   gclean plot          plot/*
%   gclean distribmodel  distribmodel
%   gclean gibbs         model_gibbs
%
%  This function is meant to avoid accidental deleting of files
%  reads .last_dir_out in the working directory
%
%  Falk Amelung, June 2007
%  modified 7/2008 so that all followers are removed by default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check remove string
%

if ~exist('.last_dir_out','file') error('.last_dir_out not found -- user error -- exiting'); end
f1=fopen('.last_dir_out','r'); last_dir_out=fscanf(f1,'%s'); fclose(f1);
[pathstr,fname,ext]  = fileparts(last_dir_out) ;   
if isempty(fname)      error('no name found -- exiting');  end

 
keepfiles = {'model_gibbs.mat'};  %need special treatment because the filematch for 'model*.mat' will otherwise remove it

if strcmp(what,'all') what = 'basemap'; keepfiles={}; end

% First string is the identifier (on command line),  second is filematch
list      = cell(2,12);                              i=1;
list(:,i) = {'basemap'};                             i=i+1;
list(:,i) = {'igram'};                               i=i+1;
list(:,i) = {'motion'};                              i=i+1;
list(:,i) = {'timeseries'};                          i=i+1;
list(:,i) = {'Qt'            'motion*Qt'};           i=i+1;
list(:,i) = {'Grid'          'motion*Grid'};         i=i+1;
list(:,i) = {'Profile'       'motion_*_Profile'};    i=i+1;
list(:,i) = {'threed'        'motion_threedfield'};  i=i+1;
list(:,i) = {'GPS'           '*GPS*'};               i=i+1;
list(:,i) = {'model'};                               i=i+1;
list(:,i) = {'distribmodel'};                        i=i+1;
list(:,i) = {'enu'};                                 i=i+1;

i_start   = strmatch(what,list(1,:),'exact');

if isempty(i_start) && isempty(strmatch(what,{'plot' 'model_gibbs' 'gibbs'})) errordlg(sprintf('string %s is not supported',what));error('-exiting'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find matching files and remove them
if exist('arg2') && strcmp('only',arg2)
    list(:,i_start+1:end)=[];
end


 for i=i_start:length(list)
     filematch = fullfile(pathstr,fname,[list{2,i} '*.mat']);
     s         = dir(filematch);
     i         = strmatch(keepfiles,{s.name}); s(i)=[];
     for j=1:length(s)
         file  = fullfile(pathstr,fname,s(j).name);
         delete(file);
         disp(['removing... ' file]) ;
     end
 end

%
% remove plot directory for gclean plot
%
 if strcmp('plot',what) 
   targetdir=[last_dir_out filesep 'plot']; 
   if exist(targetdir,'dir')  
      disp(['removing... ' targetdir]) ;
      rmdir(targetdir,'s')
   end
 end
%
% remove model_gibbs  for gclean model_gibbs
%
 if strmatch(what,{'model_gibbs' 'gibbs'},'exact') 
    file = [last_dir_out filesep 'model_gibbs.mat'];
    if exist(file,'file')
       disp(['removing... ' file]) ;
       delete(file);
    end
 end
%
% remove directory for gclean all
%
 if strcmp('basemap',what) 
   targetdir=[last_dir_out filesep 'plot']; 
   if exist(targetdir,'dir')  
      disp(['removing... ' targetdir]) ;
      rmdir(targetdir,'s')
   end
   if exist(last_dir_out,'dir') 
       disp(['removing... ' last_dir_out]) ;
       rmdir(last_dir_out,'s');           
   end
 end
