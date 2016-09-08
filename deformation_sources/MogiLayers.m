function U=MogiLayers(m,xloc,h,mu,lam,scaleN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function U=MogiLayers(m,xloc,nu,mu)
% 
% displacments at free surface due to inflationary point source
% in an elastic medium containing N layers overlying a halfspace
%
%Inputs:
%   m = 1x4 volume source geometry (length; length; length; length^3)
%        (x-coord, y-coord, depth(+), volume change)
%   xloc = 2xs matrix of observation coordinates (length)
%   h = 1xN vector of depths to bottom of N layers
%   mu = 1x(N+1) vector of shear moduli -- (last entry is shear modulus of halfspace)
%	 lam = 1x(N+1) vector of lame constant lamda -- (last entry is lamda of halfspace)
% Note: lam = -2*nu*mu/(2*nu-1), where nu is poisson's ratio
% scaleN = a positive number (set to 1 if not included as an argument)
%          scales the number of terms in the Hankel transforms
%          i.e., if there are nominally N terms, then there will be
%          round(scaleN*N) terms. The default number of terms was chosen to give
%          accurate results for observations within about 200 km of a source
%          buried less than 40 km. 
%          Run get_scaleN_mogi.m to determine number of terms needed for a particular problem.
%
% !!!Numerical Limitations: (1) large values for h, say >50, may cause the solution to blow up. 
%            If you need layers deeper than about 50, then you will have to scale down the 
%            geometry of your problem: i.e., scale all lengths by some arbitrary value
%              (2) the solution is not valid for source depths less than 0.3
%
%Outputs:
%  U = 3xs matrix of displacements
%      U(1,:) = Ux, U(2,:) = Uy, U(3,:) = Uz 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kaj Johnson, June 2002, last modified 1/10/05
if nargin==5
    scaleN=1;
end

zs=m(3);				%depth of source (km)
D=m(3);
dV=m(4);				%volume change

%find layer containing point source
temp=zs>h;
zs_layer=sum(temp)+1;
t=zs_layer;
H=h(end);				%depth to top of half_space
NL=length(h);					%number of layers above halfspace
nu=0.5*lam./(lam+mu);

Mo=2*(lam(t)+mu(t))*(1-nu(t))*dV;  %strength of source (moment) -- Mxx=Myy=Mzz=Mo
Mxx=Mo;Myy=Mo;Mzz=Mo;

muh=mu(end);					%moduli of halfspace
lamh=lam(end);

g=lam+2*mu;
gh=lamh+2*muh;

X=xloc(1,:)-m(1);		%shift coordinates
Y=xloc(2,:)-m(2);
r=sqrt(X.^2+Y.^2);
theta=atan2(Y,X);

if zs<1
	N=round(1200*scaleN);
	kmax=30;
else
   N=round(400*scaleN);
   kmax=9;
end

k=linspace(0.0001,kmax,N);

UR=zeros(1,N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate Bessel coefficients

for j=1:N
   A4x4=zeros(4,4);
   P4x4=zeros(4,4);
   
for q=1:NL
      
A4x4=[0 k(j) 1/mu(q) 0;-k(j)*lam(q)/g(q) 0 0 1/g(q);4*k(j)^2*mu(q)*(lam(q)+mu(q))/g(q) 0 0 k(j)*lam(q)/g(q);0 0 -k(j) 0];

%propagator matrix
if q==1
   z=0;
else
   z=h(q-1);
end
z0=h(q);
C3=-(sinh(k(j)*(z-z0))-k(j)*(z-z0)*cosh(k(j)*(z-z0)))/(2*k(j)^3);
C2=k(j)*(z-z0)*sinh(k(j)*(z-z0))/(2*k(j)^2);
C1=(3*sinh(k(j)*(z-z0))-k(j)*(z-z0)*cosh(k(j)*(z-z0)))/(2*k(j));
C0=(2*cosh(k(j)*(z-z0))-k(j)*(z-z0)*sinh(k(j)*(z-z0)))/2;
P4x4(:,:,q)=C3*A4x4^3+C2*A4x4^2+C1*A4x4+C0*eye(4);
end %q

%propagator from source to top of layer containing source
if zs_layer==1
   z=0;
else
   z=h(zs_layer-1);
end
z0=D;
C3=-(sinh(k(j)*(z-z0))-k(j)*(z-z0)*cosh(k(j)*(z-z0)))/(2*k(j)^3);
C2=k(j)*(z-z0)*sinh(k(j)*(z-z0))/(2*k(j)^2);
C1=(3*sinh(k(j)*(z-z0))-k(j)*(z-z0)*cosh(k(j)*(z-z0)))/(2*k(j));
C0=(2*cosh(k(j)*(z-z0))-k(j)*(z-z0)*sinh(k(j)*(z-z0)))/2;
A4x4=[0 k(j) 1/mu(t) 0;-k(j)*lam(t)/g(t) 0 0 1/g(t);4*k(j)^2*mu(t)*(lam(t)+mu(t))/g(t) 0 0 k(j)*lam(t)/g(t);0 0 -k(j) 0];
Pzs=C3*A4x4^3+C2*A4x4^2+C1*A4x4+C0*eye(4);


%product of propagator matrices
sourceP=eye(4);
if zs_layer > 1
   for q=1:zs_layer-1
      sourceP=sourceP*P4x4(:,:,q);
   end
   sourceP=sourceP*Pzs;
else
   sourceP=Pzs;
end

halfspaceP=eye(4);
for q=1:NL
   halfspaceP=halfspaceP*P4x4(:,:,q);
end
if zs > H
   halfspaceP=halfspaceP*Pzs;
end


%basis vectors for solution to homogeneous equation
%4x4
d1=[1 1 -2*muh*k(j) -2*muh*k(j)]';%*exp(-k(j)*zs);
d2=[-gh/(k(j)*(lamh+muh)) muh/(k(j)*(lamh+muh)) 2*muh 0]';%'*exp(-k(j)*zs);

%point source for explosion
%4x4
F4x4=[0;Mzz/g(t);k(j)*(-(Mxx+Myy)/2+lam(t)*Mzz/g(t));0];

%calculate constants using traction free boundary condition
B=sourceP*F4x4;
Pd1=halfspaceP*d1;
Pd2=halfspaceP*d2;
MM=[Pd1(3) Pd2(3);Pd1(4) Pd2(4)];
b=B(3:4);
unknown=inv(MM)*b;
c1=unknown(1);
c2=unknown(2);

%Bessel coefficients
US(j)=-(c1*Pd1(1)+c2*Pd2(1)-B(1));
UR(j)=c1*Pd1(2)+c2*Pd2(2)-B(2);


end %j


%%%%%%%filter noise at large k
noise=UR;   
dif=sign(diff(noise));
difdif=diff(dif);
index=(abs(difdif)>0);
temp=linspace(1,N,N);
temp=temp(index);


if length(temp)>1 
   
if temp(1)>N/10   
   filt=temp(1);
else
   filt=temp(2);
end

UR(:,filt:end)=0;
US(:,filt:end)=0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inverse HT into physical space

integrand=zeros(size(r,1),size(r,2),3);
M=0;
for j=1:length(k)
%%%Uz
   %integrand=k.*UR.*exp(i*M*theta(m,n)).*besselj(M,k*r(m,n));
	Uz_integrand(:,:,j)=k(j).*UR(j).*besselj(M,k(j)*r);
  %%%Ur 
  DrJm=-k(j)*besselj(1,k(j)*r);
  Ur_integrand(:,:,j)=US(j).*DrJm;
end

Uz=-real((1/(2*pi))*trapz(k,Uz_integrand,3));

Ur(:,:,M+1)=trapz(k,Ur_integrand,3);
Ur=real((1/(2*pi))*sum(Ur,3));
Ux=Ur.*cos(theta);
Uy=Ur.*sin(theta);

U=[Ux;Uy;Uz];