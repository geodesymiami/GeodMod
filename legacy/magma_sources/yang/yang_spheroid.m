function [sph]=spheroid(a,b,c,matrl,phi,theta,P)
% Calculate spheroid parameters and save in output vector sph

lamda=matrl(1);
mu=matrl(2);
nu=matrl(3);
% Model expressions
ac=(a-c)/(a+c);
L1=log(ac);
%L1=log((a-c)/(a+c));
iia=2/a/c^2 + L1/c^3;
iiaa=2/3/a^3/c^2 + 2/a/c^4 + L1/c^5;
coef1=-2*pi*a*b^2;
Ia=coef1*iia;
Iaa=coef1*iiaa;
u=8*pi*(1-nu);
Q=3/u;
R=(1-2*nu)/u;
a11=2*R*(Ia-4*pi);
a12=-2*R*(Ia+4*pi);
a21=Q*a^2*Iaa + R*Ia - 1;
a22=-(Q*a^2*Iaa + Ia*(2*R-Q));
coef2=3*lamda+2*mu;
w=1/(a11*a22-a12*a21);
e11=(3*a22-a12)*P*w/coef2;
e22=(a11-3*a21)*P*w/coef2;
Pdila=2*mu*(e11-e22);
Pstar=lamda*e11 + 2*(lamda+mu)*e22;
a1=-2*b^2*Pdila;
b1=3*b^2*Pdila/c^2 + 2*(1-2*nu)*Pstar; % !PL version had (1-nu) in the 2nd term!

sph(1)  = a;
sph(2)  = b;
sph(3)  = c;
sph(4)  = phi;
sph(5)  = theta;
sph(6)  = Pstar;
sph(7)  = Pdila;
sph(8)  = a1;
sph(9)  = b1;
sph(10) = P;
