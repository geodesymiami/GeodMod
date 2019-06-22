%%%%%%% Script to set the paths for GeodMod %%%%%%%%%%%%%%
%% in your startup.m file you should have:
%
%disp('Setting paths for geodmod...')
%run( [ getenv('GEODMOD_HOME') filesep 'addpath_geodmod'] )
%

disp('Added to path:')

if ~contains(path,'mimiclib')
  libdir = [ getenv('GEODMOD_HOME') filesep 'mimiclib'           ];  if exist(libdir,'dir'); addpath(genpath(libdir),'-end'); logmessage(libdir); end   % geodmod/mimiclib  (in releases)
  libdir = [ fileparts(getenv('GEODMOD_HOME')) filesep 'mimiclib'];  if exist(libdir,'dir'); addpath(genpath(libdir),'-end'); logmessage(libdir); end   % mimic/mimiclib    (at RSMAS)
end

libdir = [ getenv('GEODMOD_HOME') filesep 'masterfile'          ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMOD_HOME') filesep 'preparelib'          ];  addpath(genpath(libdir),'-end'); logmessage(libdir) 
libdir = [ getenv('GEODMOD_HOME') filesep 'modellib'            ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMOD_HOME') filesep 'PlotDatalib'         ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMOD_HOME') filesep 'plotlib'             ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMOD_HOME') filesep 'maintenancelib'      ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMOD_HOME') filesep 'deformation_sources' ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMOD_HOME') filesep 'mex'                 ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMOD_HOME') filesep 'inver'               ];  addpath(genpath(libdir),'-end'); logmessage(libdir)
libdir = [ getenv('GEODMOD_HOME') filesep 'examples'            ];  addpath(genpath(libdir),'-end'); logmessage(libdir)

libdir = [ getenv('GEODMOD_HOME') filesep 'compound_dislocations']; addpath(genpath(libdir),'-end'); logmessage(libdir)

%opengl software; logmessage('opengl set to software (at RSMAS)');
CheckToolboxes();

%set(0, 'defaultfigurewindowstyle', 'docked');
