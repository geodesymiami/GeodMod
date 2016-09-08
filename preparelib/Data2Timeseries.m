function []=Data2Timeseries(opt)
%Data2Timeseries           - reads roi_pac-processed interferograms and saves igram structure
%
%  usage:  Data2Timeseries(opt)
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
        'pattern'            ,        ['PO* filesep geo*[0-9].unw'] ,      ...
        'timeseries_name'    ,        'off'                   ,            ...
        'rating'             ,        'off'                   ,            ...
        'loadopt'             ,        'off'                   ,            ...
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

dir_in=in_name ;
filematch=fullfile(dir_in,pattern);

[junk,junk,DataSet]=extract_ProjectName(dir_in);
if ~timeseries_name
    timeseries_name         = [ 'timeseries_' DataSet  '.mat'] ;
end

%
% load timeseries
%
   %if ('subset') loadopt.subset=subset;  end 
   loadopt.unw_thresh=unw_thresh; 
   loadopt.Unit=Unit; 

   [timeseries,datelist,N_timevector]=LoadData(filematch,loadopt);

% add DataSet to timeseries structure
   for i=1:length(timeseries)
       timeseries(i).DataSet  = DataSet;
   end

   %plotdataopt.CLim='Centered';
   if isfield(plotdataopt,'CLim') plotdataopt=rmfield(plotdataopt,'CLim') ; end
   plotdataopt.colorbaropt.Location=false ;
   %logplot('PlotData_wrapper',out_name,timeseries,plotdataopt);                  % plotting give error, probably because timeseries.data is single precision

%
% save data
%
  save(out_name,'timeseries')
