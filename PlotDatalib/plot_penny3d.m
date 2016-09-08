function [s]=plot_penny3D(par)
%plot_disloc3D -  plots mogi source in 3D using fill
%
%usage:  []=PlotTheModel(par,what)
% 
%        par          dislocation model parameters 
%
% Falk Amelung, July 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = [];
if isempty(par) return ; end
logmessage(sprintf('[]=%s(%s)',mfilename,inputname(1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radius  = par(4);
[x,y,z] = sphere;
x       = x.*radius ;
y       = y.*radius;
z       = z.*radius/10;

x       = par(1) + x;
y       = par(2) + y;
z       = par(3) + z;

surface(x,y,-z,'FaceColor','red'); 

s       = par(5);

hold on





