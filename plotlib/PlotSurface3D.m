function []=PlotSurface3D(dataset,modelopt,plotdataopt,plotsurface3dopt)
%PlotSurface3D -  plots data and model in 3 dimensions
%
%usage:  []=PlotSurface3D(dataset,modelopt,plotdataopt,plotsurface3dopt)
%
%        plotsurface3dopt  options for 3-D plotting
%
%        DataSet        motion data to plot (EnvD2,East,Up,...)[dataset(1)]
%        z_offset       vertical offset of surface plot (km)   [5]
%        VertExaggeration (factor)                             [3]
%        viewdir        viewing direction (az,el)              [strike+20 20]
%        TopoMedfilt    window size for lowpass filter         [off]
%        DownsampleFac  0.5 reduces image sice by half         [1]
%        alpha_val      transparency                           [0.9]
%
% Falk Amelung, July 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dir_out plot_save
defaultopt_opt=struct(                                                     ...
        'DoIt'               ,        'on'                    ,            ...
        'fname'              ,        'surface3d'             ,            ...
        'z_offset'           ,         0                      ,            ...
        'VertExaggeration'   ,         3                      ,            ...
        'TopoMedfilt'        ,         'off'                  ,            ...
        'viewdir'            ,         'off'                  ,            ...
        'DataSet'            ,         'off'                  ,            ...
        'alpha_val'          ,         0.9                    ,            ...
        'Fringe'             ,         'on'                  ,            ...
        'Cmap'               ,         'dismph'                  ,            ...
        'DownsampleFac'      ,         'off'     )             ;        
    
[plotsurface3dopt] = process_defaultoptions(plotsurface3dopt,defaultopt_opt);  
opt=plotsurface3dopt;   f=fieldnames(opt);   for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]); end
if  ~DoIt            return; end
logmessage(sprintf('[]=%s(%s)',mfilename,inputname(1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_name    = fullfile(dir_out,fname); 
if  isempty(dataset) return; end
%if ~dataset(1).SAR   return; end
if ~dataset(1).SAR  dataset=plotdataopt.basemap; Fringe=false; end    % FA 11/2008 plot DEM if there is no SAR data. TODO: Plot modelled vertical but it is not saved by default

% display distribmodel if it exists
if isfield(modelopt,'distribmodelopt')
   plotdataopt.modelopt = modelopt.distribmodelopt;
else
    plotdataopt.modelopt = modelopt;
end

logmessage('3D plotting currently only for x_unit=km, need flakes for lola'); 

plotdataopt.x_unit           = 'km';

plotdataopt.PlotType         = '3D';
plotdataopt.z_offset         = z_offset;
plotdataopt.VertExaggeration = VertExaggeration;
plotdataopt.TopoMedfilt      = TopoMedfilt;
plotdataopt.DownsampleFac    = DownsampleFac;
plotdataopt.viewdir          = viewdir;
plotdataopt.alpha_val        = alpha_val;
plotdataopt.Fringe           = Fringe;
plotdataopt.Cmap             = Cmap;


if DataSet                                             % Select dataset to plot. Use DataSet (e.g. EnvD2) if given 
    i=strmatch(DataSet,{dataset.DataSet});             % if DataSet=Up pr East read threedfield file
    if isempty(i)
        i       = strmatch(DataSet,{'East' 'Up'});
        load(fullfile(dir_out,'motion_threedfield')); 
        dataset = enu;
        i       = strmatch(DataSet,{dataset.DataSet});
    end
else
    i=1;
end


if viewdir  viewdirs={viewdir}; else [cam_position]=viewdirections3d_generate(modelopt);end

for j=1:length(cam_position)
   plotdataopt.cam_position=cam_position{j};
   logplot('PlotData',out_name,dataset(i),plotdataopt);
end

return



