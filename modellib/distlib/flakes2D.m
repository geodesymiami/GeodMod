function [ffx,ffz]=flakes(dis_geom,nhel,nvel);
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

ffx = zeros(4,nf);
ffz = zeros(4,nf);
cnt = 1;

xloc = dis_geom(1,1)/2;

for k=1:nhel
    zloc = -10;
    for m=1:nvel
         len = dis_geom(1,cnt);
         wid = dis_geom(2,cnt);
         ffx(:,cnt) = [xloc-len/2  xloc+len/2   xloc+len/2       xloc-len/2]';
         ffz(:,cnt) = [zloc-wid/2  zloc-wid/2   zloc+wid/2       zloc+wid/2]';
         zloc = zloc+wid;
         cnt=cnt+1;
    end
    xloc = xloc+len;
end

ffz = ffz-max(ffz(:));


