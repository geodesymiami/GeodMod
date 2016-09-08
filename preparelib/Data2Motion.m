function []=Data2Motion(opt)

%Data2Timeseries           - reads roi_pac-processed interferograms and saves igram structure
%
%  usage:  Data2Motion(opt)
%
%          Options given as opt.dir_in, etc.
%
%          'dir_in'            Input directory (contains e.g. RsatA3)       [default 'off']
%          'pattern'           filematch                                    [default 'PO*/geo*[0-9].unw']
%          'timeseries_name'        filename for timeseries file                      [default is obtained with extract_name_fromSOdir]
%          'rating'            'on' or 'off' for use of igram_rating.log    [default 'off']
%          'subset'            subset to load in pixel                      [default '']
%          'unw_thresh'        unwrap threshold for LoadData                [default 0.8]
%          'corfile'           coherence to mask snaphu-unwrapped data      [default 'off']
%                              if 'off' the regular geo*cor is used
%          'cormask_thresh'    threshold for coherence mask                 [default 0.2]
%          'cormask_area_open' minimum connectd pixel in coherence mask     [default 10]
%          'Unit'              unit of data                                 [default 'radian']
%          'plotdataopt'       plot options                                 [default 'off' ]
%                        
%  Loads timeseries and saves as structure.
%  FA, Feb 2008. modified from Roipac2Igram
%  
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultopt=struct(                                                         ...
        'DoIt'               ,        'on'                   ,            ...
        'in_name'            ,        'off'                   ,            ...
        'out_name'           ,        'off'                   ,            ...
        'pattern'            ,        ['PO* filesep geo*[0-9].unw']     ,            ...
        'timeseries_name'    ,        'off'                   ,            ...
        'rating'             ,        'off'                   ,            ...
        'subset'             ,        'off'                   ,            ...
        'unw_thresh'         ,         0.99                   ,            ...
        'corfile'            ,        'off'                   ,            ...
        'cormask_thresh'     ,         0.2                    ,            ...
        'cormask_area_open'  ,         10                     ,            ...
        'Unit'               ,         'radian'               ,            ...
        'plotdataopt'        ,        'off'      )            ;
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],defaultopt); end ;
[opt]=process_defaultoptions(opt,defaultopt);  %display(opt)
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end
if  ~DoIt                    return; end
if  CheckInOut('',out_name)  return; end

dir_in    = in_name                    ;  filematch = fullfile(dir_in,pattern) ;

%
% load Motion
%
   if exist('subset') loadopt.subset = subset ;  end ;  loadopt.unw_thresh = unw_thresh ;  loadopt.Unit = Unit ; 

   [motion] = LoadData(filematch,loadopt) ;
   
   logplot('PlotData',out_name,motion,plotdataopt);  save(out_name,'motion')
