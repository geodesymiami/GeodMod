function []=PlotDataModelResidual_local(dataset,modelopt,inverseopt)
%PlotDataModelResidual_local -  plots data, model and residual in local coordinate system (using PlotModel)
%
%usage:  []=PlotDataModelResidual_local(dataset,modelopt,inverseopt)
%
%        dataset      data structure
%        model        model structure
%        inverseopt   only used to generate string containing inversemodelling method (should be removed)
%
% Falk Amelung, May 2007
% Major simplification 7/2008
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dir_out
defaultopt.plotmodelopt=struct(                                          ...
        'SampleMethod'      ,   'Grid'                  ,     ...
        'Profile'           ,   'off'                   ,     ...
        'modulo'            ,   'off'                   ,     ...
        'modulo_res'        ,   'off'                   ,     ...
        'Plot'              ,   'on'       )            ;
[plotmodelopt] = process_defaultoptions(inverseopt.plotmodelopt,defaultopt.plotmodelopt);  
logmessage(sprintf('[]=%s(%s)',mfilename,inputname(1)));

plotmodelopt.plotdataopt = inverseopt.plotdataopt;
plotmodelopt.RemoveRamps = inverseopt.PhaseRamp ;                                           % FA, Feb 2008: Got an error in PlotModel for ProjectProfile
plotmodelopt.out_name    = [dir_out '/DataModelResidual_local'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp( GenerateSummary(modelopt,dataset,inverseopt) );
invpar = modelpar2invpar(modelopt.par.xy,inverseopt.objfuncopt,1) ;
PlotModel(invpar,dataset,inverseopt,plotmodelopt) ;                       %  inverseopt used for 
