function [motion]=Timeseries2Motion(opt)

%
% Transforms timeseries measurement into motion (rate or displacement)
%
% [motion]=Timeseries2Motion(opt)
%
%   timeseries:  structure with InSAR or GPS or ? time series
%
% opt can be 
%
% 'opt'  : Structure with timeseries Info
%
%       - field      : name of the field(s) to find the motion for. 
%                      Ex: opt.InSAR.field={'data'};
%                      opt.InSAR.field(2)={'datafit'};
%
%       - model     : linear(default), bilinear,...
%
%       - time_range : [yearmin yearmax] can be multiple lines for multiple time windows
%
%
% N. Gourmelen, March 2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultopt=struct(                                                            ...
    'DoIt'                                ,        'on'          ,            ...
    'plotdataopt'                         ,        'off'         ,            ...
    'field'                               ,        'data'        ,            ...
    'time_range'                          ,        'off'         ,            ...
    'model'     	                      ,        'linear'      )             ;


if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],[]); end ;
[opt]=process_defaultoptions(opt,defaultopt); 
f=fieldnames(opt); for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end
if  ~DoIt                    return; end
if  CheckInOut(in_name,out_name)  return; end

% Load data 
if exist(in_name,'file')
    load(in_name,'timeseries')
elseif exist(in_name,var')
    timeseries = in_name ;
end
    
    
%if ~field  field='data';  end
%%%%%%%%%%%%%%%%%%%%%%%%
% Store the InSAR data %
%%%%%%%%%%%%%%%%%%%%%%%%

%phase_disp_factor=convert_unit('m/yr',timeseries);


for ni=1:length(timeseries)  Data(ni,:)=timeseries(ni).(field)(:);
    if isfield(timeseries(ni),'date_years')  datevector(ni,:)=timeseries(ni).date_years;  else  datevector(ni,:)=yymmdd2y(timeseries(ni).date);  end
end

if ~time_range  time_range = [datevector(1) datevector(length(timeseries))] ;  end

for ni=1:size(time_range,1)

in_indices     = find(datevector >= time_range(ni,1)  &  datevector <= time_range(ni,2)) ;

switch (model)
case('none')
                                                           ;
    tmp.data      = timeseries(in_indices(end)).(field)-timeseries(in_indices(1)).(field) ;
    
    tmp.time_range = time_range(ni,:) ;

    tmp = add_struct(tmp,rmfield(timeseries(in_indices(1)),'data')); % sbaker: changes timeseries(1) to timeseries(in_indices(1) so the structure
                                                                   %         has the data from the first time in the series
    motion(ni)           = tmp                                     ;
    motion(ni).TotalTime = time_range(2)-time_range(1)             ;
    motion(ni).StartDate = time_range(1)                           ;


case('linear')
                                                           ;
    motion_tmp = pinv([datevector(in_indices) ones(size(in_indices,1),1)]) * Data(in_indices,:)   ;
    tmp.data   = double(reshape(motion_tmp(1,:),timeseries(ni).file_length,timeseries(ni).width)) ;
    
    tmp.time_range = time_range(ni,:) ; 

    tmp = add_struct(tmp,rmfield(timeseries(in_indices(1)),'data')); % sbaker: changes timeseries(1) to timeseries(in_indices(1) so the structure 
                                                                   %         has the data from the first time in the series
    motion(ni)           = tmp                                     ;
    motion(ni).TotalTime = time_range(2)-time_range(1)             ;
    motion(ni).StartDate = time_range(1)                           ;
    motion(ni).data      = motion(ni).data *  motion(ni).TotalTime ;
   
end

end

if plotdataopt
    logplot('PlotData',out_name,motion,plotdataopt);
end
if out_name  save(out_name,'motion') ;  end

