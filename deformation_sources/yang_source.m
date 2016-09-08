function [u]=yang_source(param,coord,nu)
%   yang_source        - displacement due to a Yang source (simplified, no topo, equal Lame const)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Surface Displacements Due to a Yang source
% Falk Amelung, Mar 2002. 
% based on Yuri Fialko's fcn_yangM.m
%
% SIMPLIFIED:  Yuri allows for topo and different lame's parameter #########
%
%  Input:
%        param    8*1 parameters for yang source
%                     param(1) :  x-coordinate of center of source (km)
%                     param(2) :  y-coordinate of center of source (km)
%                     param(3) :  depth of center of source
%                     param(4) :  excess pressure, mu*10^(-5) Pa
%                     param(5) :  major axis, km
%                     param(6) :  ratio min-axis/may-axis  [0.01:0.99]
%                     param(7) :  strike (degree)
%                     param(8) :  plunge (degree)
%
%        coord   2*N  array  with x,y coordinates of datavec
%                     where the model displacement will be computed
%        
%  Output:
%        u        3*N array with displacement coordinates   [ux uy uz]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% start of Falk's changes %%%%


param(7)=deg2rad(param(7)); param(8)=deg2rad(param(8));
param(6)=param(5)*param(6);
as=param;
x=coord(1,:);
y=coord(2,:);
matrl(1)=1 ;
matrl(2)=1 ;
matrl(3)=nu ;
tp=zeros(size(x));
%%% end of Falk's changes %%%%

% Calculate range displacement 
mu=matrl(2);
nu=matrl(3);

% Store some commonly used parameters
coeffs(1)=1/(16*mu*(1-nu));
coeffs(2)=3-4*nu;
coeffs(3)=4*(1-nu)*(1-2*nu);

[ix,iy]=size(x);
U1r=zeros(ix,iy);
U2r=zeros(ix,iy);
U1=zeros(ix,iy);
U2=zeros(ix,iy);
U3=zeros(ix,iy);

% explicitly assign source parameters
%as
 xs    = as(1);     % center x
 ys    = as(2);     % center y
 z0    = as(3);     % center depth (positive)
 P     = as(4);     % excess pressure, mu*10^(-5) Pa
 a     = as(5);     % major axis, km
 b     = as(6);     % minor axis, km
 phi   = as(7);     % strike, rad  (0-pi)
 theta = as(8);     % plunge, rad  (0-2*pi)
 xn=x-xs;
 yn=y-ys;
 e_theta(1)=sin(theta);
 e_theta(2)=cos(theta);
 cosp=cos(phi);
 sinp=sin(phi);
 c=sqrt(a^2-b^2);
 minx=min(min(x));
 maxx=max(max(x));
 miny=min(min(y));
 maxy=max(max(y));
% if xs < minx | xs > maxx | ys < miny | ys > maxy % source is outside the grid
%  P=0;
% end

% Spheroid quantities
 [sph]=yang_spheroid(a,b,c,matrl,phi,theta,P);
 %[sph]=spheroid(a,b,c,matrl,phi,theta,P); %Anieri 4/14/15

% Rotate points
 xp=xn*cosp + yn*sinp;
 yp=yn*cosp - xn*sinp;

% Calculate model at integration limits
 xi=c;
 [Up1,Up2,Up3]=yang(sph,xi,z0,xp,yp,0,matrl,e_theta,coeffs,tp);
 xi=-xi;
 [Um1,Um2,Um3]=yang(sph,xi,z0,xp,yp,0,matrl,e_theta,coeffs,tp);

% Sum
 U1r=-Up1+Um1;
 U2r=-Up2+Um2;
% Rotate horiz. displacements back to the orig. coordinate system:
 U1=U1r*cosp-U2r*sinp+U1;
 U2=U1r*sinp+U2r*cosp+U2;
 U3=Up3-Um3+U3;

%%% next line falk's change
u=[U1;U2;U3];

