function []=PlotDataModelResidual_lola(dataset,modelopt,inverseopt)
%PlotDataModelResidual_lola -  plots data, model and residual in local coordinate system (using PlotModel)
%
%usage:  []=PlotDataModelResidual_lola(dataset,model,plotdataopt)
%
%        dataset      data structure
%        model        model structure
%        inverseopt   is only used for some calculations in PlotModel. Should be eliminated.
%
% Falk Amelung, May 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dir_out
defaultopt=struct(                                          ...
    'DoIt'               ,    'on'                   ,     ...
    'SaveModelEnu'       ,    'off'     )            ;
defaultopt.plotmodelopt=struct(                                ...
    'SampleMethod'       ,   'Grid'                  ,     ...
    'modulo'             ,   'off'                   ,     ...
    'modulo_res'        ,   'off'                   ,     ...
    'Plot'                ,   'on'       )            ;

[inverseopt] = process_defaultoptions(inverseopt,defaultopt);  % This actually should ve a different opt name
logmessage(sprintf('[]=%s(%s)',mfilename,inputname(1)));
f=fieldnames(inverseopt);  for i=1:length(f) eval([char(f{i}) '= inverseopt.(f{i}) ;' ]); end
f=fieldnames(plotdataopt); for i=1:length(f) eval([char(f{i}) '= plotdataopt.(f{i});' ]); end

if strcmp(dataset(1).CoordSystem,'ProjectProfile') return ; end

invpar=modelpar2invpar(modelopt.par.xy,objfuncopt,1) ;

plotmodelopt.plotdataopt         = struct;
plotmodelopt.Plot                = 0.5; % get the modified data according to PhaseRamp,FactorLin and predictions

[data_mod,pred,enu,dataset]=PlotModel(invpar,dataset,inverseopt,plotmodelopt) ;

N_SAR = sum([dataset(:).SAR]);
N_GPS = sum([dataset(:).GPS]);

if N_GPS
    [GPSdata,GPSpred]        = datasetstruct2GPSstruct(dataset,basemap);
    plotdataopt.GPSdata = GPSdata;
    plotdataopt.GPSpred = GPSpred;
end

if  N_SAR==0
    data=plotdataopt.basemap;
    data.SAR     = false;
    if N_GPS == 2                  % GPS horz only
        data(1).GPSonly = true ;
    end
    if N_GPS == 3                  % GPS horz and vert
        data(1).GPSonly = true ;
        %data(2).GPSonly = true ;    % FA 12/08 even if there is GPSvert there should be only one plot
    end
    N_GPSonly       = sum([data(:).GPSonly]);
else
    N_GPSonly       = 0 ;
    data(1).GPSonly = false ;
end

% fill the coordinates, Unit fields

for i=1:N_SAR*5
    data(i).x_first=dataset(1).x_first;
    data(i).y_first=dataset(1).y_first;
    data(i).x_step =dataset(1).x_step ;
    data(i).y_step =dataset(1).y_step ;
    data(i).x_unit =dataset(1).x_unit ;
    data(i).Unit   =dataset(1).Unit   ;
    data(i).SAR    =dataset(1).SAR   ;
end

% fill the data fields with the data, model, residual, ramp and original data

for i=1:N_SAR
    data(i).data            = data_mod(i).data;
    data(i).DataSet         = strcat(dataset(i).DataSet,' Data') ;
    
    data(i+N_SAR).data      = pred(i).data;
    data(i+N_SAR).DataSet   = strcat(dataset(i).DataSet,' Model') ;
    
    data(i+N_SAR*2).data    = data_mod(i).data-pred(i).data;
    data(i+N_SAR*2).DataSet = strcat(dataset(i).DataSet,' Residual') ;
    
    data(i+N_SAR*3).data    = data_mod(i).data-dataset(i).data;
    data(i+N_SAR*3).DataSet = strcat(dataset(i).DataSet,' Ramp') ;
    
    data(i+N_SAR*4).data    = dataset(i).data;
    data(i+N_SAR*4).DataSet = strcat(dataset(i).DataSet,' OrigData') ;
end
%plotdataopt.colorbaropt.Title='Data Model Residual';

out_name2 = fullfile(dir_out,'DataModelResidual_lola'    ) ;
out_name4 = fullfile(dir_out,'DataModelResidualRampOrigdata_lola') ;
data_name = fullfile(dir_out,'DataModelResidual_lola_files');

% save GPS data by Bhuvan, 2019-04
if exist('GPSdata','var')
    save(data_name,'data','GPSdata','GPSpred')
else
    save(data_name,'data') %,'GPSdata','GPSpred')
end
logplot('PlotData_DataModelResidual',out_name2,data(1:N_SAR*3+N_GPSonly),plotdataopt,3)
logplot('PlotData_DataModelResidual',out_name4,data,plotdataopt,5)
