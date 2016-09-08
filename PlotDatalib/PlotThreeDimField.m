function []=PlotData_DataModelResidual(enu,opt)
% PlotData_DataModelResidual  -  calls PlotData such that Data, Model and Residual is plotted in 3 rows             
%
%  usage:  PlotData_DataModelResidual(data,opt)
%  TODO: Needs 4th row to plot profile
if ~exist('opt','var') opt=[]; end


   opt.Quakes  = 'off' ;
   opt.gpsfile = 'off' ;
   opt.Coord   = 'off';

   PlotData(enu(3),opt)
keyboard

