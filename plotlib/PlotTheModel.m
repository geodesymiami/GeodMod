function []=PlotTheModel(dataset,modelopt,opt)
%PlotTheModel -  plots east,north,up and up with vectors of model prediction
%
%usage:  []=PlotTheModel(dataset,model,plotdataopt)
%
%        model        model structuredd
%        inverseopt   is only used for some calculations in PlotModel. Should be eliminated.
%        dataset      data structure
%
% Falk Amelung, May 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dir_out
defaultopt_opt=struct(                                      ...
        'DoIt'               ,        'on'                    ,            ...
        'QuiverGridRowsCols' ,        [15 15]                 ,            ...
        'VertExaggeration'   ,         3                      ,            ...
        'SaveModelEnu'       ,        'off'      )            ;
[opt] = process_defaultoptions(opt,defaultopt_opt);  
f=fieldnames(opt);   for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]); end
if ~DoIt  return; end
logmessage(sprintf('[]=%s(%s)',mfilename,inputname(1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% First plot East,North,Up field
%
  [enu,coord,u]        = ForwardModel_forBasemap(dataset,modelopt,plotdataopt.basemap,[],plotdataopt.x_unit);

  enu_save             = enu;
  plotdataopt.modelopt = modelopt;
  out_name             = fullfile(dir_out,'ModelENU') ;
   
  logplot('PlotData_wrapper',out_name,enu,plotdataopt)

%
%  plot vertical and horizontal vectors
%
   [junk,coord,u]  = ForwardModel_forBasemap(dataset,modelopt,plotdataopt.basemap,QuiverGridRowsCols,plotdataopt.x_unit);
   vec_east        = u(1:3:length(u));
   vec_north       = u(2:3:length(u));
   
   out_name        = fullfile(dir_out,'ModelVertandVectors') ;
   plotdataopt.VectorsBlack=[coord(1,:)',coord(2,:)',vec_east,vec_north];
   plotdataopt.colorbaropt.Location='OutsideLowerRight';

   logplot('PlotData',out_name,enu(3),plotdataopt)

%
%  plot 1D cross section
%   
      
   if isstruct(plotdataopt.Profile)
     tmpplotdataopt                  = plotdataopt;
     tmpplotdataopt.PlotType         = '1D';
     tmpplotdataopt.x_unit           = 'km';
     tmpplotdataopt.VertExaggeration = VertExaggeration;
     out_name                = fullfile(dir_out,'ModelUp_1D' );
   
     if exist([dir_out filesep 'motion_threedfield.mat'],'file') load([dir_out filesep 'motion_threedfield.mat']); end

     logplot('PlotData',out_name,enu(3),tmpplotdataopt);
   end
%
% Save modelled data as an ENU roi_pac file to aid in the roi_pac simulation process
%
if SaveModelEnu
   out_name= fullfile(dir_out,'modelenu') ;
   save (out_name,'enu_save');
   
   save_enu_roi_pac(out_name, enu_save )
end
