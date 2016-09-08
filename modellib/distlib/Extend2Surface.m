function [dis_geom_new,Shift_dist3D] = Extend2Surface(dis_geom,wid_opt)
%   Extend2Surface - Extends a dislocation to the surface
%
% Input:
% dis_geom     - (10x1) dislocation [len,width,depth,dip,strik,locE,locN,0,0,1];
% wid_opt      - (1x1) Option for the fault width, 
%                 0 - extends the fault to surface (fault width increased);
%                 1 - moves fault to surface (in fault dip direction, DEFAULT);
%
% Output:
% dis_geom_new - (10x1) new dislocation 


if nargin < 2; wid_opt=1; end
if nargin ==0 | nargin > 2; help SurfProj; end

% Input variables:
dlg = dis_geom;
if size(dlg,1)~=10; dlg=dlg';end

len = dlg(1);
wid = dlg(2);
dep = dlg(3);
dip = dlg(4);
str = dlg(5);
dE  = dlg(6);
dN  = dlg(7);

% only for negative fault dips
if dip>0 | dip < -180; disp('Use negative fault dip [-180 to 0]'); end

% Horizontal Shift
Horiz_shift  = dep / tan( dip/180*pi );
Shift_dist3D = abs ( dep / sin( dip/180*pi ) );

% Rotate shift
p     = [0 Horiz_shift]';
a     = (360-str+90)/180*pi;
R     = [cos(a) -sin(a); sin(a) cos(a)];
shift = R*p;

%new dislocation
if wid_opt == 0
   wid_new = wid + Shift_dist3D;
   dis_geom_new = [len wid_new 0 dip str dE+shift(1) dN+shift(2) dlg(8:10)']';
else
   dis_geom_new = [len wid     0 dip str dE+shift(1) dN+shift(2) dlg(8:10)']';
end
