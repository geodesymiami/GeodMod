function [fx,fy,fz]=flakes(dis_geom);
%   flakes         - Produces input for fill3 to see faultslip in 3D
%
% Program to calculate 3D cornerpoints of a rectangular fault so
% it may be plotted and visualized in 3D.  This can be particularly
% useful when distributed slip is estimated.  Use FILL3 to display
%
% Input:		dis_geom		 - (7xn) or (10xn) matrix with geometries of n faults
% output:	[fx,fy,xz]   - Input for program FILL3.  The three matrices of size
%									(4xn) have the cornerpoints of the faultpatches.
%
% e.g. dis_geom = [len,width,depth,dip,strik,locE,locN,0,0,1]';
% SJ Oct 8 2000

nf = size(dis_geom,2);

fx = zeros(4,nf);
fy = zeros(4,nf);
fz = zeros(4,nf);

for l = 1:nf
  len = dis_geom(1,l); 
  wid = dis_geom(2,l); 
  dep = dis_geom(3,l);
  dip = (dis_geom(4,l))/180*pi;
  str = dis_geom(5,l); 
  dx  = dis_geom(6,l); 
  dy  = dis_geom(7,l);

  p = [-len/2  len/2        	len/2       -len/2;
    		0      0 	     wid*cos(dip)  wid*cos(dip);
    		dep     dep			  dep-wid*sin(dip)  dep-wid*sin(dip)];

  a = (360-str+90)/180*pi;

  R = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
  k = R*p;
  kk = k + [dx;dy;0]*ones(1,4);
  
  fx(:,l) = kk(1,:)';
  fy(:,l) = kk(2,:)';
  fz(:,l) = kk(3,:)';
end
  

