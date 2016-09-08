function [pm]=patchfault(m,i,j);
%   patchfault     - discretizes a fault model into i*j distinct patches.
%PATCHFAULT    [pm]=patchfault(m,i,j)
%
%
%INPUTS:
%
%    m = 1x7 vector defined as follows
%
%       m(1) = fault length along the strike direction (km)
%       m(2) = fault width in dip direction (km)
%       m(3) = depth of lower edge of fault (km)
%       m(4) = dip angle, from the horizontal (degrees) 
%       m(5) = strike, clockwise from N (degrees)
%       m(6) = East offset of midpoint of lower edge from origin (km)
%       m(7) = North offset of midpoint of lower edge from origin (km)
%
%    i = number of patches along fault length
%
%    j = number of patches along fault width
%
%OUTPUTS:
%
%    pm = mx7 matrix of patch models, where m=i*j
%
%Examples:
% m = [100, 10, 10, 90, 0, 0, 0]
% i = 2
% j = 2
% 
% pm(1,:) = [50, 5, 5, 90, 0, 0, -25]
% pm(2,:) = [50, 5, 10, 90, 0, 0,-25]
% pm(3,:) = [50, 5, 5, 90, 0, 0, 25]
% pm(4,:) = [50, 5, 10, 90, 0, 0, 25]
%
% m = [100, 10, 10, 45, 75, 10, 10]
% i = 2
% j = 2
% 
% pm(1,:) = [50, 5, 6.4645, 45, 75, -15.0632, 6.9446]
% pm(2,:) = [50, 5, 10, 45, 75, -14.1481, 3.5295]
% pm(3,:) = [50, 5, 6.4645, 45, 75, 33.2331, 19.8855]
% pm(4,:) = [50, 5, 10, 45, 75, 34.1481, 16.4705]
 
%Set constants

	dip=m(4)*pi/180;
	strike=-m(5)*pi/180;
	sin_dip=sin(dip);
	cos_dip=cos(dip);
	iw=m(1)/i;
	jw=m(2)/j;
	is=(1:i);
	js=(1:j)';
	n=i*j;
	c1=-m(2)*cos_dip;
	c2=0.5*(m(1)+iw);
	c3=m(3)-j*jw*sin_dip;
	
%Calculate midpoints, depths of each patch
	
	p=[cos_dip*(jw*js-m(2))]*ones(1,i);
        q=ones(j,1)*[(iw*is)-0.5*(m(1)+iw)];
	r=[m(3)-jw*sin_dip*(j-js)]*ones(1,i);
	mp=[p(:),q(:),r(:)];

%Adjust midpoints for strike
	
	R=[cos(strike),-sin(strike),0;sin(strike),cos(strike),0;0,0,1];
	mp=mp*R';

%Adjust midpoints for offset from origin

	mp(:,1)=mp(:,1)+m(6);
	mp(:,2)=mp(:,2)+m(7);

%Form patch-models

	pm(:,1)=ones(n,1)*iw;
	pm(:,2)=ones(n,1)*jw;
	pm(:,3)=mp(:,3);
	pm(:,4:5)=ones(n,1)*m(4:5);
	pm(:,6:7)=mp(:,1:2);

