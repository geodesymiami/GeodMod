function [insarlist]=MakeSAR(opt)
%MakeSAR           - reads roi_pac-processed interferograms and saves igram structure
%
%  usage:  MakeSAR(opt)
%
%          Options given as opt.dir_in, etc.
%
%          'dir_in'            Input directory (contains e.g. RsatA3)       [default 'off']
%          'subset'            subset to load in pixel                      [default '']
%
%  Part of the geodmod suite
%  V1.0  Falk Amelung, June 2007
%
%   Modifications: N. Gourmelen - Feb 2008: Added testing of input type based on unit and direct reading of motion input
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dir_out;
defaultopt=struct(                                                   ...
        'DoIt'                       ,    'on'      ,                ...
        'subset'                     ,        'off'                   ,            ...
        'pattern'                    ,    ['PO*' filesep 'geo*[0-9].unw']  ,     ...
        'plotdataopt'                ,    'off'     )                ;
defaultopt.loadopt=struct(                                           ...
        'origin'                     ,    'roi_pac'  )               ;
defaultopt.samplesaropt=struct(                                      ...
         'plotdataopt'                ,   'off'                 ,     ...
         'manualSampled'              ,    {{''}}  )                  ;

defaultopt.data2igramopt.plotdataopt        =    [] ;
defaultopt.igram2motionopt.plotdataopt        =    [] ;
defaultopt.timeseries2motionopt.plotdataopt    =    [] ;
defaultopt.samplesaropt.plotdataopt           =    [] ;
defaultopt.motion2threedfieldopt.plotdataopt  =    [] ;
defaultopt.dir_inlist                         =    {} ;  % For some reason the assignment of {} is
defaultopt.motion_list                        =    {} ;  % otherwise not possible (matlab bug?)
defaultopt.qt_list                            =    {} ;
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],defaultopt); end ;
[opt]=process_defaultoptions(opt,defaultopt);  
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end
if  ~DoIt || ~length(dir_inlist) insarlist={};      return; end
%if  CheckInOut('',out_name)  return; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data2igramopt.plotdataopt         = process_defaultoptions(data2igramopt.plotdataopt        ,plotdataopt);
igram2motionopt.plotdataopt       = process_defaultoptions(igram2motionopt.plotdataopt      ,plotdataopt);
timeseries2motionopt.plotdataopt  = process_defaultoptions(timeseries2motionopt.plotdataopt ,plotdataopt);
motion2threedfieldopt.plotdataopt = process_defaultoptions(motion2threedfieldopt.plotdataopt,plotdataopt);
samplesaropt.plotdataopt.basemap  = plotdataopt.basemap;
samplesaropt.plotdataopt.CLim     = plotdataopt.CLim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Make the motion files for each dataset         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(dir_inlist)

    [igram_list{i},timeseries_list{i},motion_list{i},qt_list{i}] = GenerateNames(dir_inlist{i},dir_out,samplesaropt);
    
    loadopt.subset                =  subset             ;
    data2igramopt.in_name         =  dir_inlist{i}      ;
    data2igramopt.out_name        =  igram_list{i}      ;
    data2igramopt.pattern         =  pattern            ;
    data2igramopt.loadopt         =  loadopt            ;

    data2timeseriesopt.in_name    =  dir_inlist{i}      ;
    data2timeseriesopt.out_name   =  timeseries_list{i} ;
    data2timeseriesopt.pattern    =  pattern            ;
    data2timeseriesopt.loadopt    =  loadopt            ;

    igram2motionopt.in_name       =  igram_list{i}      ;
    igram2motionopt.out_name      =  motion_list{i}     ;

    timeseries2motionopt.in_name  =  timeseries_list{i} ;
    timeseries2motionopt.out_name =  motion_list{i}     ;

    loadopt.subset                =  subset             ;
    data2motionopt.in_name        =  dir_inlist{i}      ;
    data2motionopt.out_name       =  motion_list{i}     ;
    data2motionopt.pattern        =  pattern            ;
    data2motionopt.loadopt        =  loadopt            ;
    
    if iscell(data2igramopt.cormask_thresh)                                      % FA 1/2011: switching input to cell arrays
       data2igramopt.cormask_thresh_value  =  data2igramopt.cormask_thresh{i};
    end
    
    samplesaropt.in_name          = motion_list{i}     ;
    samplesaropt.out_name         = char(qt_list{i})   ;
    
    if ~isempty(samplesaropt.manualSampled{1})                                   % FA 1/2011: If manualSampled(1) ~= {''} (the default) then use cell contents
       samplesaropt.manualSampledFile= samplesaropt.manualSampled{i};          
    end

    data_type = CheckDatatype(data2igramopt);

    switch data_type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'radian'
            Data2Igram(data2igramopt)               ;
            Igram2Motion(igram2motionopt)           ;
        case 'm'
            Data2Timeseries(data2timeseriesopt)     ;
            Timeseries2Motion(timeseries2motionopt) ;
        case 'm/yr'
            Data2Motion(data2motionopt)             ; % FA 12/2008 field StartTime not yet introduced in Data2Motion as I forgot what it is for.
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SampleSar(samplesaropt)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  ~CheckInOut('',[ dir_out filesep 'plot' filesep 'TimeCoverage.pdf' ])
    logplot('PlotTimeCoverage',[dir_out filesep 'TimeCoverage'],motion_list);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now process the datasets simultaneouly        %%
motion2threedfieldopt.in_list  = motion_list;

Motion2ThreeDfield(motion2threedfieldopt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
insarlist=qt_list;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
