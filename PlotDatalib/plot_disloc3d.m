function [s]=plot_disloc3D(par,what)
%plot_disloc3D -  plots dislocation in 3D using fill
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
%logmessage(sprintf('[]=%s(%s)',mfilename,inputname(1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load slipcol    % Load colormap
viewdir = [par(5)+10 20];        % default azimuth clode to strike

% convert par into the patchmodel format if necessary

if size(par,2) == 10                             % already in N-by-10 format
       pm       = par ;
elseif mod(size(par,1),10)==0 && size(par,2)==1    % in 10*N-by-1 format
       N_disloc = size(par,1)/10;   
       pm       = reshape(par',10,N_disloc)';
end
    
% 
switch what
    case 'ss'
        s       = pm(:,8);
        partype = 'Strike-Slip';
    case 'ds'
        s       = pm(:,9);
        partype = 'Dip-Slip';
    case 'op'
        s       = pm(:,10);
        partype = 'Opening'; 
    case 'mag'
        s       = sqrt(pm(:,8).^2 + pm(:,9).^2);
        rake    = 180-rad2deg(atan2(pm(:,9),pm(:,8)));
        partype = 'Slip Magnitude';
    otherwise
        errordlg('argument needs to be ss,ds,op,mag'); error('user error -- exiting');     
end
    
[fx,fy,fz] = flakes(pm');
fill3(fx,fy,-fz,s','LineStyle','-'); 
%axis image;

str = sprintf( 'Predicted %s Distribution',partype); 
%str= sprintf( 'Predicted %s Distribution, roughness = %4.2f cm/km',ParType{i},roughness*100); 

title(str); 

hold on





