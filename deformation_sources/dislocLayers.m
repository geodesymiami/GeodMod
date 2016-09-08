function U=dislocLayers(m,xloc,d,mu,lam,scaleN)

%U=LayeredGreens(m,d,mu,lam,xloc)
% Calculates the surface displacements for a dislocation in a layered
% elastic halfspace using propagator matrix methods.
%%INPUTS:
%fault model (standard Okada paramterization -- must be column vector):
% m(1) = length (km) (do not give distances in meters)
% m(2) = width (km)
% m(3) = depth to down-dip edge (km)F
% m(4) = dip (degrees)
% m(5) = strike (degrees)
% m(6) = east position of down-dip edge (km)
% m(7) = north position (km)
% m(8) = strike-slip (m) (units can be anythiing, actually -- units of
%        output displacements will be same as slip input units since output
%        displacements scale linearly with slip)
% m(9) = dip-slip (m)
% m(10) = tensile (m)
% d = 1x(M-1) vector of depths to bottom of layer -- M is number of layers
%     (including halfspace); note: bottom layer is halfspace, thus no depth is specified
% mu,lam = 1xM vector of normalized Lame constants (normalized with value of mu in top layer)
%          mu=lam for poisson ratio = 0.25, M is number of layers (including halfspace)
%          It is IMPORTANT that the Lame constants used are relative to (normalized by)
%          the value in the top layer
% Example: An earth model with two layers overlying a halfspace, the depths
% of the two layers are 10 and 20 km with mu=3*10^10 in first layer and
% mu=5*10^10 in second layer and mu=7*10^10 below the second layer in the
% halfspace with poisson ratio in all layers equal to 0.25. Set d=[10 20],
% mu=[1 5/3 7/3], lam=[1 5/3, 7/3]
%
% xloc = 2xn matrix of n station coordinates (km) -- first row x, second row y
% scaleN = a positive number (set to 1 if not included as an argument)
%          scales the number of terms in the Hankel transforms
%          i.e., if there are nominally N terms, then there will be
%          round(scaleN*N) terms. Run get_scaleN.m to determine number of
%          terms needed for a particular problem.
%%OUTPUTS:
%U is a 3xn matrix of displacements (first row - east, second row - north, third row - up) 
%
%Kaj M. Johnson, Stanford and UC Berkeley, last updated 1/10/05

if nargin==5
    scaleN=1;  %set scaleN to default
end

xy=xloc'; 
%initialize basis vectors to zero
G1=zeros(3*size(xy,1),1);
G2=zeros(3*size(xy,1),1);
G3=zeros(3*size(xy,1),1);

%define components of slip
ss=-m(8);
ds=m(9);
ten=m(10);


%convert parameters to be consistent with Okada
if m(4)<=90 & m(4)>=-90
   dipdir = (m(5) + 90)*pi/180;
else
   dipdir = (m(5) - 90)*pi/180;
end

offset = abs(m(2).*cos(m(4)*pi/180));
EastOff = offset.*sin(dipdir);
NorthOff = offset.*cos(dipdir);
xy(:,1)=xy(:,1)-m(6)+EastOff;
xy(:,2)=xy(:,2)-m(7)+NorthOff;

L=m(1);
W=m(2);
dip=m(4);
strike=m(5);
D=m(3)-W*sin(dip*pi/180);

%shift and rotate coordinates such that fault has zero strike and centered at origin
R=[cos(strike*pi/180) -sin(strike*pi/180);sin(strike*pi/180) cos(strike*pi/180)];	%rotation matrix about z-axis
rotcoords=R*[xy(:,1)';xy(:,2)'];
xy(:,1)=rotcoords(1,:)';
xy(:,2)=rotcoords(2,:)';
strike=0;

Xcoord=xy(:,2);  %switch x and y
Ycoord=xy(:,1);

muh=mu(end);
lamh=lam(end);
for k=1:length(d)
   g(k)=lam(k)+2*mu(k);
end
gh=lamh+2*muh;
g=[g gh];

%position of point sources
%setup spacing of point sources for Gaussian quadrature
%number of point sources near ground surface
   

NL=ceil(m(1));
NW=ceil(m(2));

u=(1:NL-1)./sqrt((2*(1:NL-1)).^2-1);
	[vc,bp]=eig(diag(u,-1)+diag(u,1));
	[bp,k]=sort(diag(bp));
	a=-L/2;b=L/2;
	wL=(2*vc(1,k).^2)*(b-a)/2;
	xp=(a+b)/2+(b-a)/2*bp;

	u=(1:NW-1)./sqrt((2*(1:NW-1)).^2-1);
	[vc,bp]=eig(diag(u,-1)+diag(u,1));
	[bp,k]=sort(diag(bp));
   
   a=0;b=W;


   wW=(2*vc(1,k).^2)*(b-a)/2;
	yp=(a+b)/2+(b-a)/2*bp;
   
[xp,yp]=meshgrid(xp,yp);   
  
xp=xp(:)';
yp=yp(:)';
zp=zeros(size(xp));


wW=repmat(wW,size(xy,1),1);
wW=repmat(wW,[1 1 NL]);
wL=reshape(wL,[1 1 NL]);
wL=repmat(wL,size(xy,1),NW);

%mesh points
%rotate into true position
%rotate plane to true dip
R=[1 0 0;0 cos(dip*pi/180) -sin(dip*pi/180);0 sin(dip*pi/180) cos(dip*pi/180)];	%rotation matrix x-axis
Xp=R*[xp;yp;zp];
%rotate to true strike
R=[cos(strike*pi/180) -sin(strike*pi/180) 0;sin(strike*pi/180) cos(strike*pi/180) 0;0 0 1];	%rotation matrix about z-axis
Xp=R*Xp;
xpos=Xp(1,:);
ypos=Xp(2,:);
zpos=Xp(3,:);
xpos=reshape(xpos,NW,NL);
ypos=reshape(ypos,NW,NL);
zpos=zpos(1:NW);


for psW=1:NW

      
zs=D+zpos(psW);		%depth of point source


ZS=zs;
if zs<0
   disp('Warning: fault extends above ground surface')
end



if zs>=.5 
   sourceloop=1;
else 
   sourceloop=2;
end



%initialize variables
for jj=1:3
ur_0_1{jj}=[];ur_0_2{jj}=[];ur_0_3{jj}=[];
us_0_1{jj}=[];us_0_2{jj}=[];us_0_3{jj}=[];

ur_n1_1{jj}=[];ur_n1_2{jj}=[];ur_n1_3{jj}=[];
us_n1_1{jj}=[];us_n1_2{jj}=[];us_n1_3{jj}=[];
ut_n1_1{jj}=[];ut_n1_2{jj}=[];ut_n1_3{jj}=[];

ur_p1_1{jj}=[];ur_p1_2{jj}=[];ur_p1_3{jj}=[];
us_p1_1{jj}=[];us_p1_2{jj}=[];us_p1_3{jj}=[];
ut_p1_1{jj}=[];ut_p1_2{jj}=[];ut_p1_3{jj}=[];

ur_n2_1{jj}=[];ur_n2_2{jj}=[];ur_n2_3{jj}=[];
us_n2_1{jj}=[];us_n2_2{jj}=[];us_n2_3{jj}=[];
ut_n2_1{jj}=[];ut_n2_2{jj}=[];ut_n2_3{jj}=[];

ur_p2_1{jj}=[];ur_p2_2{jj}=[];ur_p2_3{jj}=[];
us_p2_1{jj}=[];us_p2_2{jj}=[];us_p2_3{jj}=[];
ut_p2_1{jj}=[];ut_p2_2{jj}=[];ut_p2_3{jj}=[];
end
UR_0_1=[];UR_0_2=[];UR_0_3=[];
US_0_1=[];US_0_2=[];US_0_3=[];

UR_n1_1=[];UR_n1_2=[];UR_n1_3=[];
US_n1_1=[];US_n1_2=[];US_n1_3=[];
UT_n1_1=[];UT_n1_2=[];UT_n1_3=[];

UR_p1_1=[];UR_p1_2=[];UR_p1_3=[];
US_p1_1=[];US_p1_2=[];US_p1_3=[];
UT_p1_1=[];UT_p1_2=[];UT_p1_3=[];

UR_n2_1=[];UR_n2_2=[];UR_n2_3=[];
US_n2_1=[];US_n2_2=[];US_n2_3=[];
UT_n2_1=[];UT_n2_2=[];UT_n2_3=[];

UR_p2_1=[];UR_p2_2=[];UR_p2_3=[];
US_p2_1=[];US_p2_2=[];US_p2_3=[];
UT_p2_1=[];UT_p2_2=[];UT_p2_3=[];

for J=1:sourceloop
   
if ZS>=7   
   
Filter=1;   
NN=round(300*scaleN);
N=round(300*scaleN);

kdepth=[7 11 15 25 30];
kmaxdepth=[3 2 1 .7 .5];
index=sum(ZS>=kdepth);
kmax=kmaxdepth(index);

k=linspace(0.001,kmax,NN);
K=k;

elseif ZS>=4   
   
Filter=1;   
NN=round(400*scaleN);
N=round(400*scaleN);

kdepth=[4 7 11 15 25 30];
kmaxdepth=[4 3 2 1 .7 .5];
index=sum(ZS>=kdepth);
kmax=kmaxdepth(index);

k=linspace(0.001,kmax,NN);
K=k;

elseif ZS>=2  
   
Filter=1;   
NN=round(600*scaleN);
N=round(600*scaleN);

kdepth=[2 3 4 7 11 15 25 30];
kmaxdepth=[6 5 4 3 2 1 .7 .5];
index=sum(ZS>=kdepth);
kmax=kmaxdepth(index);

k=linspace(0.001,kmax,NN);
K=k;

elseif ZS>=.5   
   
Filter=1;   
NN=round(800*scaleN);
N=round(800*scaleN);

kdepth=[.5 1 2 3 4 7 11 15 25 30];
kmaxdepth=[20 8 6 5 4 3 2 1 .7 .5];
index=sum(ZS>=kdepth);
kmax=kmaxdepth(index);

k=linspace(0.001,kmax,NN);
K=k;


elseif ZS>=.12
   
Filter=0;   
%NN=7200;
NN=round(1200*scaleN);

kdepth=[.1 .2 .3];
kmaxdepth=[100 60 50];
index=sum(ZS>=kdepth);
kmax=kmaxdepth(index);

K=linspace(.001,kmax,NN);
if J==1
   index=logical(K<=30);
   N=sum(index);
   k=K(index);
else
   index=logical(K>30);
   N=sum(index);
   k=K(index)/10;
   zs=10*zs;
end

else
    
Filter=0;   
NN=round(4000*scaleN);

kdepth=[0 .07 .12];
kmaxdepth=[180 150 100];
index=sum(ZS>=kdepth);
kmax=kmaxdepth(index);

K=linspace(.001,kmax,NN);
if J==1
   index=logical(K<=30);
   N=sum(index);
   k=K(index);
else
   index=logical(K>30);
   N=sum(index);
   k=K(index)/10;
   zs=10*zs;
end

end %if zs


kstore{psW}=K;
Nstore{psW}=NN;

%find layer containing point source
temp=zs>d;
zs_layer=sum(temp)+1;
t=zs_layer;

[M1,M2,M3]=momtensor_inverse(strike,dip,lam(t),mu(t));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate Bessel coefficients
UR_0_1=zeros(1,N);US_0_1=zeros(1,N);
UR_n1_1=zeros(1,N);US_n1_1=zeros(1,N);UT_n1_1=zeros(1,N);
UR_p1_1=zeros(1,N);US_p1_1=zeros(1,N);UT_p1_1=zeros(1,N);
UR_n2_1=zeros(1,N);US_n2_1=zeros(1,N);UT_n2_1=zeros(1,N);
UR_p2_1=zeros(1,N);US_p2_1=zeros(1,N);UT_p2_1=zeros(1,N);

UR_0_2=zeros(1,N);US_0_2=zeros(1,N);
UR_n1_2=zeros(1,N);US_n1_2=zeros(1,N);UT_n1_2=zeros(1,N);
UR_p1_2=zeros(1,N);US_p1_2=zeros(1,N);UT_p1_2=zeros(1,N);
UR_n2_2=zeros(1,N);US_n2_2=zeros(1,N);UT_n2_2=zeros(1,N);
UR_p2_2=zeros(1,N);US_p2_2=zeros(1,N);UT_p2_2=zeros(1,N);

UR_0_3=zeros(1,N);US_0_3=zeros(1,N);
UR_n1_3=zeros(1,N);US_n1_3=zeros(1,N);UT_n1_3=zeros(1,N);
UR_p1_3=zeros(1,N);US_p1_3=zeros(1,N);UT_p1_3=zeros(1,N);
UR_n2_3=zeros(1,N);US_n2_3=zeros(1,N);UT_n2_3=zeros(1,N);
UR_p2_3=zeros(1,N);US_p2_3=zeros(1,N);UT_p2_1=zeros(1,N);
bignum=0;
j=0;
while bignum ==0 
    
    j=j+1;
   
    
   A4x4=zeros(4,4);
   P4x4=zeros(4,4,length(d));
   
[sourceP4x4,sourceP2x2,halfspaceP4x4,halfspaceP2x2]=getprop(d,k,j,mu,lam,g,zs,zs_layer);
   
%basis vectors for solution to homogeneous equation
H=d(end);
if zs>H;
   z0=zs;
else
   z0=H;
end

   a2=[-.5/muh -.5/muh k(j) k(j)]';
   b2=[1/(2*k(j)*(lamh+muh));-(lamh+2*muh)/(2*muh*k(j)*(lamh+muh));0;1];
   d1=a2;
   d2=b2+a2*z0;
   
   
%2x2
h=[1/muh -abs(k(j))]';




%calculate constants using traction free boundary condition
Pd1=halfspaceP4x4*d1;
Pd2=halfspaceP4x4*d2;
MM=[Pd1(3) Pd2(3);Pd1(4) Pd2(4)];
Ph=halfspaceP2x2*h;




%% make source vectors and get bessel coefficients



if ss~=0
   [bignum,UR_0_1(:,j),US_0_1(:,j), ...
         UR_n1_1(:,j),US_n1_1(:,j),UT_n1_1(:,j), ...
         UR_p1_1(:,j),US_p1_1(:,j),UT_p1_1(:,j), ...
         UR_n2_1(:,j),US_n2_1(:,j),UT_n2_1(:,j), ...
         UR_p2_1(:,j),US_p2_1(:,j),UT_p2_1(:,j)]=getbesselco(bignum,MM,Pd1,Pd2,Ph,sourceP4x4,...
   sourceP2x2,M1,g,t,k,j,lam,mu);
end


if ds~=0
[bignum,UR_0_2(:,j),US_0_2(:,j), ...
         UR_n1_2(:,j),US_n1_2(:,j),UT_n1_2(:,j), ...
         UR_p1_2(:,j),US_p1_2(:,j),UT_p1_2(:,j), ...
         UR_n2_2(:,j),US_n2_2(:,j),UT_n2_2(:,j), ...
         UR_p2_2(:,j),US_p2_2(:,j),UT_p2_2(:,j)]=getbesselco(bignum,MM,Pd1,Pd2,Ph,sourceP4x4,sourceP2x2,M2,g,t,k,j,lam,mu);
end



if ten~=0
[bignum,UR_0_3(:,j),US_0_3(:,j), ...
         UR_n1_3(:,j),US_n1_3(:,j),UT_n1_3(:,j), ...
         UR_p1_3(:,j),US_p1_3(:,j),UT_p1_3(:,j), ...
         UR_n2_3(:,j),US_n2_3(:,j),UT_n2_3(:,j), ...
         UR_p2_3(:,j),US_p2_3(:,j),UT_p2_3(:,j)]=getbesselco(bignum,MM,Pd1,Pd2,Ph,sourceP4x4,...
   sourceP2x2,M3,g,t,k,j,lam,mu);
end

if zs<1; bignum=0; end
if j==N; bignum=1;end

end %while (j)




%%%%%%%filter noise at large k
if ds==0
	noise=imag(UR_n2_1(1,:));   
else
   noise=imag(UR_n1_2(1,:));
end

dif=diff(noise,10);
maxval=max(abs(noise(1:50)));
%index=dif>maxval*10^-3;
index=dif>.2*10^-3;
index(1:50)=0;

dummy=linspace(1,length(k),length(k));
filt=dummy(index);

if ~isempty(filt)

if ss~=0
UR_0_1(:,filt:end)=0;
US_0_1(:,filt:end)=0;

UR_n1_1(:,filt:end)=0;
US_n1_1(:,filt:end)=0;
UT_n1_1(:,filt:end)=0;

UR_p1_1(:,filt:end)=0;
US_p1_1(:,filt:end)=0;
UT_p1_1(:,filt:end)=0;

UR_n2_1(:,filt:end)=0;
US_n2_1(:,filt:end)=0;
UT_n2_1(:,filt:end)=0;

UR_p2_1(:,filt:end)=0;
US_p2_1(:,filt:end)=0;
UT_p2_1(:,filt:end)=0;
end

if ds~=0
UR_0_2(:,filt:end)=0;
US_0_2(:,filt:end)=0;

UR_n1_2(:,filt:end)=0;
US_n1_2(:,filt:end)=0;
UT_n1_2(:,filt:end)=0;

UR_p1_2(:,filt:end)=0;
US_p1_2(:,filt:end)=0;
UT_p1_2(:,filt:end)=0;

UR_n2_2(:,filt:end)=0;
US_n2_2(:,filt:end)=0;
UT_n2_2(:,filt:end)=0;

UR_p2_2(:,filt:end)=0;
US_p2_2(:,filt:end)=0;
UT_p2_2(:,filt:end)=0;
end

if ten~=0
UR_0_3(:,filt:end)=0;
US_0_3(:,filt:end)=0;

UR_n1_3(:,filt:end)=0;
US_n1_3(:,filt:end)=0;
UT_n1_3(:,filt:end)=0;

UR_p1_3(:,filt:end)=0;
US_p1_3(:,filt:end)=0;
UT_p1_3(:,filt:end)=0;

UR_n2_3(:,filt:end)=0;
US_n2_3(:,filt:end)=0;
UT_n2_3(:,filt:end)=0;

UR_p2_3(:,filt:end)=0;
US_p2_3(:,filt:end)=0;
UT_p2_3(:,filt:end)=0;
end

end

if ss~=0
ur_0_1{J}=UR_0_1;
us_0_1{J}=US_0_1;

ur_n1_1{J}=UR_n1_1;
us_n1_1{J}=US_n1_1;
ut_n1_1{J}=UT_n1_1;

ur_p1_1{J}=UR_p1_1;
us_p1_1{J}=US_p1_1;
ut_p1_1{J}=UT_p1_1;

ur_n2_1{J}=UR_n2_1;
us_n2_1{J}=US_n2_1;
ut_n2_1{J}=UT_n2_1;

ur_p2_1{J}=UR_p2_1;
us_p2_1{J}=US_p2_1;
ut_p2_1{J}=UT_p2_1;
end

if ds~=0
ur_0_2{J}=UR_0_2;
us_0_2{J}=US_0_2;

ur_n1_2{J}=UR_n1_2;
us_n1_2{J}=US_n1_2;
ut_n1_2{J}=UT_n1_2;

ur_p1_2{J}=UR_p1_2;
us_p1_2{J}=US_p1_2;
ut_p1_2{J}=UT_p1_2;

ur_n2_2{J}=UR_n2_2;
us_n2_2{J}=US_n2_2;
ut_n2_2{J}=UT_n2_2;

ur_p2_2{J}=UR_p2_2;
us_p2_2{J}=US_p2_2;
ut_p2_2{J}=UT_p2_2;
end

if ten~=0
ur_0_3{J}=UR_0_3;
us_0_3{J}=US_0_3;

ur_n1_3{J}=UR_n1_3;
us_n1_3{J}=US_n1_3;
ut_n1_3{J}=UT_n1_3;

ur_p1_3{J}=UR_p1_3;
us_p1_3{J}=US_p1_3;
ut_p1_3{J}=UT_p1_3;

ur_n2_3{J}=UR_n2_3;
us_n2_3{J}=US_n2_3;
ut_n2_3{J}=UT_n2_3;

ur_p2_3{J}=UR_p2_3;
us_p2_3{J}=US_p2_3;
ut_p2_3{J}=UT_p2_3;
end


end %J



if ss~=0
	   UR_0_1s{psW}=[ur_0_1{1} ur_0_1{2} ur_0_1{3}];
      US_0_1s{psW}=[us_0_1{1} us_0_1{2} us_0_1{3}];
      
  	   UR_n1_1s{psW}=[ur_n1_1{1} ur_n1_1{2} ur_n1_1{3}];
      US_n1_1s{psW}=[us_n1_1{1} us_n1_1{2} us_n1_1{3}];
      UT_n1_1s{psW}=[ut_n1_1{1} ut_n1_1{2} ut_n1_1{3}];
      
  	   UR_p1_1s{psW}=[ur_p1_1{1} ur_p1_1{2} ur_p1_1{3}];
      US_p1_1s{psW}=[us_p1_1{1} us_p1_1{2} us_p1_1{3}];
      UT_p1_1s{psW}=[ut_p1_1{1} ut_p1_1{2} ut_p1_1{3}];
      
  	   UR_n2_1s{psW}=[ur_n2_1{1} ur_n2_1{2} ur_n2_1{3}];
      US_n2_1s{psW}=[us_n2_1{1} us_n2_1{2} us_n2_1{3}];
      UT_n2_1s{psW}=[ut_n2_1{1} ut_n2_1{2} ut_n2_1{3}];
      
  	   UR_p2_1s{psW}=[ur_p2_1{1} ur_p2_1{2} ur_p2_1{3}];
      US_p2_1s{psW}=[us_p2_1{1} us_p2_1{2} us_p2_1{3}];
      UT_p2_1s{psW}=[ut_p2_1{1} ut_p2_1{2} ut_p2_1{3}];
end
   
if ds~=0
	   UR_0_2s{psW}=[ur_0_2{1} ur_0_2{2} ur_0_2{3}];
      US_0_2s{psW}=[us_0_2{1} us_0_2{2} us_0_2{3}];
      
  	   UR_n1_2s{psW}=[ur_n1_2{1} ur_n1_2{2} ur_n1_2{3}];
      US_n1_2s{psW}=[us_n1_2{1} us_n1_2{2} us_n1_2{3}];
      UT_n1_2s{psW}=[ut_n1_2{1} ut_n1_2{2} ut_n1_2{3}];
      
  	   UR_p1_2s{psW}=[ur_p1_2{1} ur_p1_2{2} ur_p1_2{3}];
      US_p1_2s{psW}=[us_p1_2{1} us_p1_2{2} us_p1_2{3}];
      UT_p1_2s{psW}=[ut_p1_2{1} ut_p1_2{2} ut_p1_2{3}];
      
  	   UR_n2_2s{psW}=[ur_n2_2{1} ur_n2_2{2} ur_n2_2{3}];
      US_n2_2s{psW}=[us_n2_2{1} us_n2_2{2} us_n2_2{3}];
      UT_n2_2s{psW}=[ut_n2_2{1} ut_n2_2{2} ut_n2_2{3}];
      
  	   UR_p2_2s{psW}=[ur_p2_2{1} ur_p2_2{2} ur_p2_2{3}];
      US_p2_2s{psW}=[us_p2_2{1} us_p2_2{2} us_p2_2{3}];
      UT_p2_2s{psW}=[ut_p2_2{1} ut_p2_2{2} ut_p2_2{3}];
end

if ten~=0
	   UR_0_3s{psW}=[ur_0_3{1} ur_0_3{2} ur_0_3{3}];
      US_0_3s{psW}=[us_0_3{1} us_0_3{2} us_0_3{3}];
      
  	   UR_n1_3s{psW}=[ur_n1_3{1} ur_n1_3{2} ur_n1_3{3}];
      US_n1_3s{psW}=[us_n1_3{1} us_n1_3{2} us_n1_3{3}];
      UT_n1_3s{psW}=[ut_n1_3{1} ut_n1_3{2} ut_n1_3{3}];
      
  	   UR_p1_3s{psW}=[ur_p1_3{1} ur_p1_3{2} ur_p1_3{3}];
      US_p1_3s{psW}=[us_p1_3{1} us_p1_3{2} us_p1_3{3}];
      UT_p1_3s{psW}=[ut_p1_3{1} ut_p1_3{2} ut_p1_3{3}];
      
  	   UR_n2_3s{psW}=[ur_n2_3{1} ur_n2_3{2} ur_n2_3{3}];
      US_n2_3s{psW}=[us_n2_3{1} us_n2_3{2} us_n2_3{3}];
      UT_n2_3s{psW}=[ut_n2_3{1} ut_n2_3{2} ut_n2_3{3}];
      
  	   UR_p2_3s{psW}=[ur_p2_3{1} ur_p2_3{2} ur_p2_3{3}];
      US_p2_3s{psW}=[us_p2_3{1} us_p2_3{2} us_p2_3{3}];
      UT_p2_3s{psW}=[ut_p2_3{1} ut_p2_3{2} ut_p2_3{3}];
end

end %psW


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inverse HT into physical space


for psW=1:NW

    zs=D+zpos(psW);
    
   for psL=1:NL


xshift=repmat(xpos(psW,psL),length(Xcoord),1);
yshift=repmat(ypos(psW,psL),length(Xcoord),1);

X=Xcoord-xshift;
Y=Ycoord-yshift;

r=sqrt(X.^2+Y.^2);
theta=atan2(Y,X);

k=kstore{psW};
N=Nstore{psW};     

rr=r;
kk=k;
r=repmat(r,1,length(k));
k=repmat(k,length(rr),1);
kr=k.*r;
theta_orig=theta;
theta=repmat(theta,1,N);

if ss~=0
UR_0_1=repmat(UR_0_1s{psW},length(rr),1);
US_0_1=repmat(US_0_1s{psW},length(rr),1);

UR_n1_1=repmat(UR_n1_1s{psW},length(rr),1);
US_n1_1=repmat(US_n1_1s{psW},length(rr),1);
UT_n1_1=repmat(UT_n1_1s{psW},length(rr),1);

UR_p1_1=repmat(UR_p1_1s{psW},length(rr),1);
US_p1_1=repmat(US_p1_1s{psW},length(rr),1);
UT_p1_1=repmat(UT_p1_1s{psW},length(rr),1);

UR_n2_1=repmat(UR_n2_1s{psW},length(rr),1);
US_n2_1=repmat(US_n2_1s{psW},length(rr),1);
UT_n2_1=repmat(UT_n2_1s{psW},length(rr),1);

UR_p2_1=repmat(UR_p2_1s{psW},length(rr),1);
US_p2_1=repmat(US_p2_1s{psW},length(rr),1);
UT_p2_1=repmat(UT_p2_1s{psW},length(rr),1);
end

if ds~=0
UR_0_2=repmat(UR_0_2s{psW},length(rr),1);
US_0_2=repmat(US_0_2s{psW},length(rr),1);
UR_n1_2=repmat(UR_n1_2s{psW},length(rr),1);
US_n1_2=repmat(US_n1_2s{psW},length(rr),1);
UT_n1_2=repmat(UT_n1_2s{psW},length(rr),1);

UR_p1_2=repmat(UR_p1_2s{psW},length(rr),1);
US_p1_2=repmat(US_p1_2s{psW},length(rr),1);
UT_p1_2=repmat(UT_p1_2s{psW},length(rr),1);

UR_n2_2=repmat(UR_n2_2s{psW},length(rr),1);
US_n2_2=repmat(US_n2_2s{psW},length(rr),1);
UT_n2_2=repmat(UT_n2_2s{psW},length(rr),1);

UR_p2_2=repmat(UR_p2_2s{psW},length(rr),1);
US_p2_2=repmat(US_p2_2s{psW},length(rr),1);
UT_p2_2=repmat(UT_p2_2s{psW},length(rr),1);
end

if ten~=0
UR_0_3=repmat(UR_0_3s{psW},length(rr),1);
US_0_3=repmat(US_0_3s{psW},length(rr),1);
UR_n1_3=repmat(UR_n1_3s{psW},length(rr),1);
US_n1_3=repmat(US_n1_3s{psW},length(rr),1);
UT_n1_3=repmat(UT_n1_3s{psW},length(rr),1);

UR_p1_3=repmat(UR_p1_3s{psW},length(rr),1);
US_p1_3=repmat(US_p1_3s{psW},length(rr),1);
UT_p1_3=repmat(UT_p1_3s{psW},length(rr),1);

UR_n2_3=repmat(UR_n2_3s{psW},length(rr),1);
US_n2_3=repmat(US_n2_3s{psW},length(rr),1);
UT_n2_3=repmat(UT_n2_3s{psW},length(rr),1);

UR_p2_3=repmat(UR_p2_3s{psW},length(rr),1);
US_p2_3=repmat(US_p2_3s{psW},length(rr),1);
UT_p2_3=repmat(UT_p2_3s{psW},length(rr),1);
end


M=0;

b0=besselj(0,kr);
b1=besselj(1,kr);
b2=-b0+2./kr.*b1;
b3=-b1+4./kr.*b2;

DrJm=-k.*b1;

   if ss~=0
	Uz_integrand_1=k.*UR_0_1.*b0;
	Ur_integrand_1=US_0_1.*DrJm;
	end
   
   if ds~=0
   Uz_integrand_2=k.*UR_0_2.*b0;
   Ur_integrand_2=US_0_2.*DrJm;
	end
   
   if ten~=0
   Uz_integrand_3=k.*UR_0_3.*b0;
   Ur_integrand_3=US_0_3.*DrJm;
   end


    

if ss~=0
Uz_1(:,1)=(1/(2*pi))*trapz(kk,Uz_integrand_1,2);
Ur_1(:,1)=real((1/(2*pi))*trapz(kk,Ur_integrand_1,2));
Utheta_1(:,1)=zeros(length(rr),1);
end

if ds~=0
Uz_2(:,1)=(1/(2*pi))*trapz(kk,Uz_integrand_2,2);
Ur_2(:,1)=real((1/(2*pi))*trapz(kk,Ur_integrand_2,2));
Utheta_2(:,1)=zeros(length(rr),1);
end

if ten~=0
Uz_3(:,1)=(1/(2*pi))*trapz(kk,Uz_integrand_3,2);
Ur_3(:,1)=real((1/(2*pi))*trapz(kk,Ur_integrand_3,2));
Utheta_3(:,1)=zeros(length(rr),1);
end

M=-1;
DrJm=.5*k.*(b2-b0);
   
   if ss~=0
	Uz_integrand_1=-k.*UR_n1_1.*b1.*exp(i*M*theta);
	Ur_integrand_1=(US_n1_1.*DrJm - i*M*UT_n1_1.*(1./r).*b1).*exp(i*M*theta);
	end
   
   if ds~=0
	Uz_integrand_2=-k.*UR_n1_2.*b1.*exp(i*M*theta);
	Ur_integrand_2=(US_n1_2.*DrJm - i*M*UT_n1_2.*(1./r).*b1).*exp(i*M*theta);
	end
   
   if ten~=0
 	Uz_integrand_3=-k.*UR_n1_3.*b1.*exp(i*M*theta);
	Ur_integrand_3=(US_n1_3.*DrJm - i*M*UT_n1_3.*(1./r).*b1).*exp(i*M*theta);
	end


if ss~=0
Uz_1(:,2)=(1/(2*pi))*trapz(kk,Uz_integrand_1,2);
Ur_1(:,2)=real((1/(2*pi))*trapz(kk,Ur_integrand_1,2));
end

if ds~=0
Uz_2(:,2)=(1/(2*pi))*trapz(kk,Uz_integrand_2,2);
Ur_2(:,2)=real((1/(2*pi))*trapz(kk,Ur_integrand_2,2));
Utheta_2(:,2)=zeros(length(rr),1);
end

if ten~=0
Uz_3(:,2)=(1/(2*pi))*trapz(kk,Uz_integrand_3,2);
Ur_3(:,2)=real((1/(2*pi))*trapz(kk,Ur_integrand_3,2));
Utheta_3(:,2)=zeros(length(rr),1);
end

M=1;
DrJm=.5*k.*(b0-b2);

   if ss~=0
	Uz_integrand_1=k.*UR_p1_1.*b1.*exp(i*M*theta);
	Ur_integrand_1=(US_p1_1.*DrJm + i*M*UT_p1_1.*(1./r).*b1).*exp(i*M*theta);
	end
   
   if ds~=0
	Uz_integrand_2=k.*UR_p1_2.*b1.*exp(i*M*theta);
	Ur_integrand_2=(US_p1_2.*DrJm + i*M*UT_p1_2.*(1./r).*b1).*exp(i*M*theta);
	end
   
   if ten~=0
	Uz_integrand_3=k.*UR_p1_3.*b1.*exp(i*M*theta);
	Ur_integrand_3=(US_p1_3.*DrJm + i*M*UT_p1_3.*(1./r).*b1).*exp(i*M*theta);
	end


if ss~=0
Uz_1(:,3)=(1/(2*pi))*trapz(kk,Uz_integrand_1,2);
Ur_1(:,3)=real((1/(2*pi))*trapz(kk,Ur_integrand_1,2));
Utheta_1(:,3)=zeros(length(rr),1);
end

if ds~=0
Uz_2(:,3)=(1/(2*pi))*trapz(kk,Uz_integrand_2,2);
Ur_2(:,3)=real((1/(2*pi))*trapz(kk,Ur_integrand_2,2));
Utheta_2(:,3)=zeros(length(rr),1);
end

if ten~=0
Uz_3(:,3)=(1/(2*pi))*trapz(kk,Uz_integrand_3,2);
Ur_3(:,3)=real((1/(2*pi))*trapz(kk,Ur_integrand_3,2));
Utheta_3(:,3)=zeros(length(rr),1);
end


M=-2;
DrJm=.5*k.*(b1-b3);
   
   if ss~=0
	Uz_integrand_1=k.*UR_n2_1.*b2.*exp(i*M*theta);
	Ur_integrand_1=(US_n2_1.*DrJm + i*M*UT_n2_1.*(1./r).*b2).*exp(i*M*theta);
	Utheta_integrand_1=(i*M*US_n2_1.*(1./r).*b2 - UT_n2_1.*DrJm).*exp(i*M*theta);
	end
   
   if ds~=0
	Uz_integrand_2=k.*UR_n2_2.*b2.*exp(i*M*theta);
	Ur_integrand_2=(US_n2_2.*DrJm + i*M*UT_n2_2.*(1./r).*b2).*exp(i*M*theta);
	Utheta_integrand_2=(i*M*US_n2_2.*(1./r).*b2 - UT_n2_2.*DrJm).*exp(i*M*theta);
  	end  
  
   if ten~=0
	Uz_integrand_3=k.*UR_n2_3.*b2.*exp(i*M*theta);
	Ur_integrand_3=(US_n2_3.*DrJm + i*M*UT_n2_3.*(1./r).*b2).*exp(i*M*theta);
	Utheta_integrand_3=(i*M*US_n2_3.*(1./r).*b2 - UT_n2_3.*DrJm).*exp(i*M*theta);
	end  

if ss~=0
Uz_1(:,4)=(1/(2*pi))*trapz(kk,Uz_integrand_1,2);
Ur_1(:,4)=real((1/(2*pi))*trapz(kk,Ur_integrand_1,2));
Utheta_1(:,4)=real((1/(2*pi))*trapz(kk,Utheta_integrand_1,2));
end

if ds~=0
Uz_2(:,4)=(1/(2*pi))*trapz(kk,Uz_integrand_2,2);
Ur_2(:,4)=real((1/(2*pi))*trapz(kk,Ur_integrand_2,2));
Utheta_2(:,4)=real((1/(2*pi))*trapz(kk,Utheta_integrand_2,2));
end

if ten~=0
Uz_3(:,4)=(1/(2*pi))*trapz(kk,Uz_integrand_3,2);
Ur_3(:,4)=real((1/(2*pi))*trapz(kk,Ur_integrand_3,2));
Utheta_4(:,3)=real((1/(2*pi))*trapz(kk,Utheta_integrand_3,2));
end


M=2;
DrJm=.5*k.*(b1-b3);
   
   if ss~=0
	Uz_integrand_1=k.*UR_p2_1.*b2.*exp(i*M*theta);
	Ur_integrand_1=(US_p2_1.*DrJm + i*M*UT_p2_1.*(1./r).*b2).*exp(i*M*theta);
	Utheta_integrand_1=(i*M*US_p2_1.*(1./r).*b2 - UT_p2_1.*DrJm).*exp(i*M*theta);
	end
   
   if ds~=0
	Uz_integrand_2=k.*UR_p2_2.*b2.*exp(i*M*theta);
	Ur_integrand_2=(US_p2_2.*DrJm + i*M*UT_p2_2.*(1./r).*b2).*exp(i*M*theta);
	Utheta_integrand_2=(i*M*US_p2_2.*(1./r).*b2 - UT_p2_2.*DrJm).*exp(i*M*theta);
  	end  
  
   if ten~=0
	Uz_integrand_3=k.*UR_p2_3.*b2.*exp(i*M*theta);
	Ur_integrand_3=(US_p2_3.*DrJm + i*M*UT_p2_3.*(1./r).*b2).*exp(i*M*theta);
	Utheta_integrand_3=(i*M*US_p2_3.*(1./r).*b2 - UT_p2_3.*DrJm).*exp(i*M*theta);
	end  

if ss~=0
Uz_1(:,5)=(1/(2*pi))*trapz(kk,Uz_integrand_1,2);
Ur_1(:,5)=real((1/(2*pi))*trapz(kk,Ur_integrand_1,2));
Utheta_1(:,5)=real((1/(2*pi))*trapz(kk,Utheta_integrand_1,2));
end

if ds~=0
Uz_2(:,5)=(1/(2*pi))*trapz(kk,Uz_integrand_2,2);
Ur_2(:,5)=real((1/(2*pi))*trapz(kk,Ur_integrand_2,2));
Utheta_2(:,5)=real((1/(2*pi))*trapz(kk,Utheta_integrand_2,2));
end

if ten~=0
Uz_3(:,5)=(1/(2*pi))*trapz(kk,Uz_integrand_3,2);
Ur_3(:,5)=real((1/(2*pi))*trapz(kk,Ur_integrand_3,2));
Utheta_5(:,3)=real((1/(2*pi))*trapz(kk,Utheta_integrand_3,2));
end

if ss~=0
Gz_1(:,psW,psL)=sum(Uz_1,2);
UR=sum(Ur_1,2);
UTHETA=sum(Utheta_1,2);
Urx=UR.*cos(theta_orig);
Ury=UR.*sin(theta_orig);
Uthetax=UTHETA.*cos(theta_orig+pi/2);
Uthetay=UTHETA.*sin(theta_orig+pi/2);
Gx_1(:,psW,psL)=Urx+Uthetax;
Gy_1(:,psW,psL)=Ury+Uthetay;

if zs<1   %remove aliasing in far-field
index=r(:,1)-m(7)>100;
Gx_1(index,psW,psL)=0;
Gy_1(index,psW,psL)=0;
Gz_1(index,psW,psL)=0;
end

end


if ds~=0
Gz_2(:,psW,psL)=sum(Uz_2,2);
UR=sum(Ur_2,2);
UTHETA=sum(Utheta_2,2);
Urx=UR.*cos(theta_orig);
Ury=UR.*sin(theta_orig);
Uthetax=UTHETA.*cos(theta_orig+pi/2);
Uthetay=UTHETA.*sin(theta_orig+pi/2);
Gx_2(:,psW,psL)=Urx+Uthetax;
Gy_2(:,psW,psL)=Ury+Uthetay;

if zs<2   %remove aliasing in far-field
index=r(:,1)-m(7)>100;
Gx_2(index,psW,psL)=0;
Gy_2(index,psW,psL)=0;
Gz_2(index,psW,psL)=0;
end

end

if ten~=0
Gz_3(:,psW,psL)=sum(Uz_3,2);
UR=sum(Ur_3,2);
UTHETA=sum(Utheta_3,2);
Urx=UR.*cos(theta_orig);
Ury=UR.*sin(theta_orig);
Uthetax=UTHETA.*cos(theta_orig+pi/2);
Uthetay=UTHETA.*sin(theta_orig+pi/2);
Gx_3(:,psW,psL)=Urx+Uthetax;
Gy_3(:,psW,psL)=Ury+Uthetay;

if zs<1.5   %remove aliasing in far-field
index=r(:,1)-m(7)>100;
Gx_3(index,psW,psL)=0;
Gy_3(index,psW,psL)=0;
Gz_3(index,psW,psL)=0;
end

end


end %psW
end %psL





if ss~=0
   
   Gx_1=wW.*wL.*Gx_1;
   Gy_1=wW.*wL.*Gy_1;
   Gz_1=wW.*wL.*Gz_1;
   
   Gx_1=sum(sum(Gx_1,2),3);
	Gy_1=sum(sum(Gy_1,2),3);
	Gz_1=sum(sum(Gz_1,2),3);

	G1=zeros(3*length(Gy_1),1);
	G1(1:3:end)=Gy_1;
	G1(2:3:end)=Gx_1;
	G1(3:3:end)=-Gz_1;
end

if ds~=0
   
   Gx_2=wW.*wL.*Gx_2;
   Gy_2=wW.*wL.*Gy_2;
   Gz_2=wW.*wL.*Gz_2;
   
   Gx_2=sum(sum(Gx_2,2),3);
	Gy_2=sum(sum(Gy_2,2),3);
	Gz_2=sum(sum(Gz_2,2),3);

	G2=zeros(3*length(Gy_2),1);
	G2(1:3:end)=Gy_2;
	G2(2:3:end)=Gx_2;
   G2(3:3:end)=-Gz_2;
end

if ten~=0
   
   Gx_3=wW.*wL.*Gx_3;
   Gy_3=wW.*wL.*Gy_3;
   Gz_3=wW.*wL.*Gz_3;
   

   Gx_3=sum(sum(Gx_3,2),3);
	Gy_3=sum(sum(Gy_3,2),3);
	Gz_3=sum(sum(Gz_3,2),3);

	G3=zeros(3*length(Gy_3),1);
	G3(1:3:end)=Gy_3;
	G3(2:3:end)=Gx_3;
   G3(3:3:end)=-Gz_3;

end


if ss==0
   G1=zeros(size(xloc,2)*3,1);
end
if ds==0
   G2=zeros(size(xloc,2)*3,1);
end
if ten==0
   G3=zeros(size(xloc,2)*3,1);
end


%change sign to be consistent with Okada's definition of slip
G2=-G2;

%%output displacments in matrix U
U=zeros(3,length(G1)/3);
%scale displacements with slip
temp=ss*G1+ds*G2+ten*G3;
U(1,:)=temp(1:3:end)';
U(2,:)=temp(2:3:end)';
U(3,:)=temp(3:3:end)';      %U is 3xn matrix of displacements


%rotate back into true position
strike=-m(5);
R=[cos(strike*pi/180) -sin(strike*pi/180);sin(strike*pi/180) cos(strike*pi/180)];	%rotation matrix about z-axis
rotdisps=R*[U(1,:);U(2,:)];
U(1,:)=rotdisps(1,:);
U(2,:)=rotdisps(2,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sub functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M1,M2,M3]=momtensor_inverse(strike,dip,lam,mu)

strike=strike*pi/180;
dip=dip*pi/180;

Vs=[cos(strike) sin(strike) 0]';
Vd=[-cos(dip)*sin(strike) cos(dip)*cos(strike) sin(dip)]';
Vnorm=cross(Vd,Vs);

%moment tensor
%strike slip component
M1(1,1)=Vs(1)*Vnorm(1)*(lam+2*mu)+Vs(2)*Vnorm(2)*lam+Vs(3)*Vnorm(3)*lam;
M1(2,2)=Vs(1)*Vnorm(1)*lam+Vs(2)*Vnorm(2)*(lam+2*mu)+Vs(3)*Vnorm(3)*lam;
M1(3,3)=Vs(1)*Vnorm(1)*lam+Vs(2)*Vnorm(2)*lam+Vs(3)*Vnorm(3)*(lam+2*mu);
M1(1,2)=Vs(1)*Vnorm(2)*mu+Vs(2)*Vnorm(1)*mu;
M1(1,3)=Vs(1)*Vnorm(3)*mu+Vs(3)*Vnorm(1)*mu;
M1(2,3)=Vs(2)*Vnorm(3)*mu+Vs(3)*Vnorm(2)*mu;
M1(2,1)=M1(1,2);M1(3,1)=M1(1,3);M1(3,2)=M1(2,3);
%dip slip component
M2(1,1)=Vd(1)*Vnorm(1)*(lam+2*mu)+Vd(2)*Vnorm(2)*lam+Vd(3)*Vnorm(3)*lam;
M2(2,2)=Vd(1)*Vnorm(1)*lam+Vd(2)*Vnorm(2)*(lam+2*mu)+Vd(3)*Vnorm(3)*lam;
M2(3,3)=Vd(1)*Vnorm(1)*lam+Vd(2)*Vnorm(2)*lam+Vd(3)*Vnorm(3)*(lam+2*mu);
M2(1,2)=Vd(1)*Vnorm(2)*mu+Vd(2)*Vnorm(1)*mu;
M2(1,3)=Vd(1)*Vnorm(3)*mu+Vd(3)*Vnorm(1)*mu;
M2(2,3)=Vd(2)*Vnorm(3)*mu+Vd(3)*Vnorm(2)*mu;
M2(2,1)=M2(1,2);M2(3,1)=M2(1,3);M2(3,2)=M2(2,3);
%tensile component
M3(1,1)=Vnorm(1)*Vnorm(1)*(lam+2*mu)+Vnorm(2)*Vnorm(2)*lam+Vnorm(3)*Vnorm(3)*lam;
M3(2,2)=Vnorm(1)*Vnorm(1)*lam+Vnorm(2)*Vnorm(2)*(lam+2*mu)+Vnorm(3)*Vnorm(3)*lam;
M3(3,3)=Vnorm(1)*Vnorm(1)*lam+Vnorm(2)*Vnorm(2)*lam+Vnorm(3)*Vnorm(3)*(lam+2*mu);
M3(1,2)=Vnorm(1)*Vnorm(2)*mu+Vnorm(2)*Vnorm(1)*mu;
M3(1,3)=Vnorm(1)*Vnorm(3)*mu+Vnorm(3)*Vnorm(1)*mu;
M3(2,3)=Vnorm(2)*Vnorm(3)*mu+Vnorm(3)*Vnorm(2)*mu;
M3(2,1)=M3(1,2);M3(3,1)=M3(1,3);M3(3,2)=M3(2,3);

temp=M1./max(max(abs(M1)));
index=logical(abs(temp)<10^-13);
M1(index)=0;

temp=M2./max(max(abs(M2)));
index=logical(abs(temp)<10^-13);
M2(index)=0;

temp=M3./max(max(abs(M3)));
index=logical(abs(temp)<10^-13);
M3(index)=0;


function [sourceP4x4,sourceP2x2,halfspaceP4x4,halfspaceP2x2]=getprop(d,k,j,mu,lam,g,zs,zs_layer)


%subroutine to generate propagator matrices
t=zs_layer;
H=d(end);
NL=length(d);
for q=1:length(d)
      
A4x4=[0 k(j) 1/mu(q) 0;-k(j)*lam(q)/g(q) 0 0 1/g(q);4*k(j)^2*mu(q)*(lam(q)+mu(q))/g(q) 0 0 k(j)*lam(q)/g(q);0 0 -k(j) 0];


if q==1
   z=0;
else
   z=d(q-1);
end
z0=d(q);
C3=-(sinh(k(j)*(z-z0))-k(j)*(z-z0)*cosh(k(j)*(z-z0)))/(2*k(j)^3);
C2=k(j)*(z-z0)*sinh(k(j)*(z-z0))/(2*k(j)^2);
C1=(3*sinh(k(j)*(z-z0))-k(j)*(z-z0)*cosh(k(j)*(z-z0)))/(2*k(j));
C0=(2*cosh(k(j)*(z-z0))-k(j)*(z-z0)*sinh(k(j)*(z-z0)))/2;


P4x4(:,:,q)=C3*A4x4^3+C2*A4x4^2+C1*A4x4+C0*eye(4);
P2x2(1,:,q)=[cosh((z-z0)*abs(k(j))) (1/(mu(q)*abs(k(j))))*sinh((z-z0)*abs(k(j)))];
P2x2(2,:,q)=[mu(q)*abs(k(j))*sinh((z-z0)*abs(k(j))) cosh((z-z0)*abs(k(j)))];
end %q



%propagator from source to top of layer containing source
if zs_layer==1
   z=0;
else
   z=d(zs_layer-1);
end
z0=zs;
C3=-(sinh(k(j)*(z-z0))-k(j)*(z-z0)*cosh(k(j)*(z-z0)))/(2*k(j)^3);
C2=k(j)*(z-z0)*sinh(k(j)*(z-z0))/(2*k(j)^2);
C1=(3*sinh(k(j)*(z-z0))-k(j)*(z-z0)*cosh(k(j)*(z-z0)))/(2*k(j));
C0=(2*cosh(k(j)*(z-z0))-k(j)*(z-z0)*sinh(k(j)*(z-z0)))/2;
A4x4=[0 k(j) 1/mu(t) 0;-k(j)*lam(t)/g(t) 0 0 1/g(t);4*k(j)^2*mu(t)*(lam(t)+mu(t))/g(t) 0 0 k(j)*lam(t)/g(t);0 0 -k(j) 0];
P4x4zs=C3*A4x4^3+C2*A4x4^2+C1*A4x4+C0*eye(4);
P2x2zs(1,:)=[cosh((z-z0)*abs(k(j))) (1/(mu(t)*abs(k(j))))*sinh((z-z0)*abs(k(j)))];
P2x2zs(2,:)=[mu(t)*abs(k(j))*sinh((z-z0)*abs(k(j))) cosh((z-z0)*abs(k(j)))];



%product of propagator matrices
halfspaceP4x4=eye(4);
halfspaceP2x2=eye(2);
for q=1:NL
   halfspaceP4x4=halfspaceP4x4*P4x4(:,:,q);
   halfspaceP2x2=halfspaceP2x2*P2x2(:,:,q);
end
if zs > H
   halfspaceP4x4=halfspaceP4x4*P4x4zs;
   halfspaceP2x2=halfspaceP2x2*P2x2zs;
end


sourceP4x4=eye(4);
sourceP2x2=eye(2);
if zs_layer > 1
   for q=1:zs_layer-1
      sourceP4x4=sourceP4x4*P4x4(:,:,q);
      sourceP2x2=sourceP2x2*P2x2(:,:,q);
   end
   sourceP4x4=sourceP4x4*P4x4zs;
   sourceP2x2=sourceP2x2*P2x2zs;
else
   sourceP4x4=P4x4zs;
   sourceP2x2=P2x2zs;
end



function [bignum,UR_0,US_0,UR_n1,US_n1,UT_n1,UR_p1,US_p1,UT_p1,...
      UR_n2,US_n2,UT_n2,UR_p2,US_p2,UT_p2]=getbesselco(bignum,MM,Pd1,Pd2,...
   Ph,sourceP4x4,sourceP2x2,M,g,t,k,j,lam,mu)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%subroutine to generate bessel coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mxx=M(1,1);
Myy=M(2,2);
Mzz=M(3,3);
Mxy=M(1,2);
Mxz=M(1,3);
Myz=M(2,3);

%point source
%order = 0
F1_0=0;
F2_0=Mzz/g(t);
F3_0=k(j)*((-Mxx-Myy)/2+lam(t)*Mzz/g(t));

%order = -1
F1_n1=-1/2*(Mxz+i*Myz)/mu(t);
F2_n1=0;
F3_n1=0;
f1_n1=1/2*(Myz-i*Mxz)/mu(t);
f2_n1=0;

%order = 1
F1_p1=-1/2*(-Mxz+i*Myz)/mu(t);
F2_p1=0;
F3_p1=0;
f1_p1=-1/2*(Myz+i*Mxz)/mu(t);
f2_p1=0;

%order = -2
F1_n2=0;
F2_n2=0;
F3_n2=k(j)*(-(Myy-Mxx)/4+i*Mxy/2);
f1_n2=0;
f2_n2=k(j)*(i*(Mxx-Myy)/4-1/2*Mxy);

%order = -2
F1_p2=0;
F2_p2=0;
F3_p2=k(j)*(-(Myy-Mxx)/4-i*Mxy/2);
f1_p2=0;
f2_p2=k(j)*(-i*(Mxx-Myy)/4-1/2*Mxy);


      %order = 0
      F4x4=[F1_0;F2_0;F3_0;0];
      B=sourceP4x4*F4x4;	
		b=B(3:4,:);
      unknown=inv(MM)*b;
		c1=unknown(1,:);
      c2=unknown(2,:);
      UR_0=(c1*Pd1(2)+c2*Pd2(2)-B(2,:))';
		US_0=(-(c1*Pd1(1)+c2*Pd2(1)-B(1,:)))';
   
      if  log10(mean(mean(MM))) > 30
          bignum=1;
      end
      
  
      
      %order = -1
      F4x4=[F1_n1;F2_n1;F3_n1;0];
      B=sourceP4x4*F4x4;	
		b=B(3:4,:);
        unknown=inv(MM)*b;
		c1=unknown(1,:);
      c2=unknown(2,:);
      UR_n1=(c1*Pd1(2)+c2*Pd2(2)-B(2,:))';
      US_n1=(-(c1*Pd1(1)+c2*Pd2(1)-B(1,:)))';
      f=[f1_n1;f2_n1]; 
		Pf=sourceP2x2*f;
		c=Pf(2,:)/Ph(2);
      UT_n1=(c*Ph(1)-Pf(1,:))';
      
      
      
      %order = 1
      F4x4=[F1_p1;F2_p1;F3_p1;0];
      B=sourceP4x4*F4x4;	
		b=B(3:4,:);
		unknown=inv(MM)*b;
		c1=unknown(1,:);
      c2=unknown(2,:);
      UR_p1=(c1*Pd1(2)+c2*Pd2(2)-B(2,:))';
      US_p1=(-(c1*Pd1(1)+c2*Pd2(1)-B(1,:)))';
      f=[f1_p1;f2_p1]; 
		Pf=sourceP2x2*f;
		c=Pf(2,:)/Ph(2);
      UT_p1=(c*Ph(1)-Pf(1,:))';
      
      %order = -2
      F4x4=[F1_n2;F2_n2;F3_n2;0];
      B=sourceP4x4*F4x4;	
		b=B(3:4,:);
		unknown=inv(MM)*b;
		c1=unknown(1,:);
      c2=unknown(2,:);
      UR_n2=(c1*Pd1(2)+c2*Pd2(2)-B(2,:))';
      US_n2=(-(c1*Pd1(1)+c2*Pd2(1)-B(1,:)))';
      f=[f1_n2;f2_n2]; 
		Pf=sourceP2x2*f;
		c=Pf(2,:)/Ph(2);
      UT_n2=(c*Ph(1)-Pf(1,:))';
      
      %order = 2
      F4x4=[F1_p2;F2_p2;F3_p2;0];
      B=sourceP4x4*F4x4;	
		b=B(3:4,:);
		unknown=inv(MM)*b;
		c1=unknown(1,:);
      c2=unknown(2,:);
      UR_p2=(c1*Pd1(2)+c2*Pd2(2)-B(2,:))';
      US_p2=(-(c1*Pd1(1)+c2*Pd2(1)-B(1,:)))';
      f=[f1_p2;f2_p2]; 
		Pf=sourceP2x2*f;
		c=Pf(2,:)/Ph(2);
      UT_p2=(c*Ph(1)-Pf(1,:))';

      

