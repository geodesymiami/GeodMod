path='.';                        % path to where the data is
gotodir=['cd ' path];
eval(gotodir);

x0    = 0;   % origin coordinates, km
y0    = 0;
dx    = 0.5;   % grid size, km
dy    = 0.5;
;
ndat=100;     % number of points in the grid for x & y
mdat=100;
xx=x0:dx:(mdat-1)*dx+x0;
yy=y0:dy:(ndat-1)*dy+y0;
[x,y]=meshgrid(xx,yy);
tp=zeros(size(x));   % no topo

as=zeros(1,8);         % model vector  (see fcn_yangM.m)
as(1)=20;
as(2)=30;
as(3)=15;
as(4)=10;
as(5)=12;
as(6)=4;
as(7)=pi/180*40;
as(8)=pi/180*30;


dhx=round(ndat/10);       % subsampled grid:
dhy=round(mdat/10);
xsub = x(1:dhx:ndat,1:dhy:mdat);
ysub = y(1:dhx:ndat,1:dhy:mdat);

% Elastic constants:
matrl(1)=1;         % 1st Lame constant
matrl(2)=matrl(1);  % shear modulus (2nd Lame constant)
matrl(3)=0.25;      % Poisson's ratio
