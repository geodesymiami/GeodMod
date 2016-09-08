function [fringe]=CalcDefaultFringe(data, PlotUnit)
%CalcDefaultFringe  - calculates fringe spacing so that 1 fringe corresponds to wavelength/2 
%
%usage: [satparameters]=extract_hardwired_satparameters(in_name,parameter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavelength = extract_hardwired_satparameters(data.DataSet,'wavelength');

switch PlotUnit
case {'m'}
     fringe = wavelength/2;
case {'m/yr'}
     fringe = wavelength/2/data.TotalTime;
case {'radian'}
     fringe = 2*pi;
otherwise
     errordlg(sprintf(['PlotUnit ' PlotUnit ' not yet supported for fringes'])) ; error('user error -- exiting');
end
