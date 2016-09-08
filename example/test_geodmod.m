function []=test_geodmod(what)
%test_geodmod    -  test geodmod distribution
%
%  usage:  test_geodmod
%  usage:  test_geodmod 'all'
%
% Function to test geodmod distribution. Assumes Data (DEM, interferograms, GPS, linefiles,
% quakefiles, etc) are available in the GEODMOD_TESTDATA directory.  Files are written to
% the GEODMOD_TESTBENCH directory and removed after completion
%
% Without any arguments it runs Darwin.min, Forward_mogi.min, MLGPSonly and
% basemap_hawaii.min
%
% For 'all' it also runs ML2002-2005, Darwin_pennyshapedcrack, Darwin_yang.
%
%  Quick testing by running:
%
%  geodmod Darwin.min
%  geodmod Forward_mogi.min
%  geodmod MLGPSonly.min
%  geodmod basemap_hawaii.min
%  geodmod ML2002-2005.min
%  geodmod Darwin_pennyshapedcrack
%  geodmod Darwin_yang
%  

  all=false ;
  if nargin==1 && (strcmp(what,'all') || strcmp(what,'All'))         % all=true if all tests are run
     all=true ;     
  end
  disp(all)


  list   = {'Darwin' 'Forward_mogi' 'MLGPSonly' 'basemap_hawaii' 'Wells'};
  others = {'ML2002-2005' 'Darwin_pennyshapedcrack' 'Darwin_yang'}; 

  if all   list={list{:} others{:}};  end

  for i=1:length(list)
      dirname  = fullfile(getenv('GEODMOD_TESTBENCH'),list{i}) ;
      filename = [getenv('GEODMODHOME') filesep 'distribution_testing' filesep list{i} '.min' ] ;
      if exist(dirname,'dir')   rmdir(dirname,'s');  end
      logmessage(['Test #' num2str(i) ', running....   geodmod ' filename])

      geodmod(filename,'plot_save=0');

      gclean all
  end

  clear all ; close all          % needed fo that no error messages in batch mode
  system('rm -f .last_dir_out');
