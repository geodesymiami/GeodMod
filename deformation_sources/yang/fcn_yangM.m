function [U1,U2,U3]=fcn_yangM(as,x,y,matrl,tp)
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
 phi   = as(7);     % plunge, rad  (0-pi)
 theta = as(8);     % strike, rad  (0-2*pi)
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

% Speroid quantities
 [sph]=spheroid(a,b,c,matrl,phi,theta,P);

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



