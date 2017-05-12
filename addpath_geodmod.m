%%%%%%% Script to set the paths for GeodMod %%%%%%%%%%%%%%
%% in your startup.m file you should have:
%
%disp('Setting paths for geodmod...')
%run( [ getenv('GEODMODHOME') filesep 'addpath_geodmod'] )
%

disp('Added to path:')

if isempty(strfind(path,'mimiclib'))
  libdir = [ getenv('GEODMODHOME') filesep 'mimiclib'           ];  if exist(libdir,'dir') addpath(genpath(libdir),'-end'); logmessage(libdir); end   % geodmod/mimiclib  (in releases)
  libdir = [ fileparts(getenv('GEODMODHOME')) filesep 'mimiclib'];  if exist(libdir,'dir') addpath(genpath(libdir),'-end'); logmessage(libdir); end   % mimic/mimiclib    (at RSMAS)
end

libdir = [ getenv('GEODMODHOME') filesep 'masterfile'          ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMODHOME') filesep 'preparelib'          ];  addpath(genpath(libdir),'-end'); logmessage(libdir) 
libdir = [ getenv('GEODMODHOME') filesep 'modellib'            ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMODHOME') filesep 'PlotDatalib'         ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMODHOME') filesep 'plotlib'             ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMODHOME') filesep 'maintenancelib'      ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMODHOME') filesep 'deformation_sources' ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMODHOME') filesep 'mex'                 ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMODHOME') filesep 'inver'               ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMODHOME') filesep 'examples'            ];  addpath(genpath(libdir),'-end'); logmessage(libdir)

%opengl software; logmessage('opengl set to software (at RSMAS)');
CheckToolboxes();

%set(0, 'defaultfigurewindowstyle', 'docked');
