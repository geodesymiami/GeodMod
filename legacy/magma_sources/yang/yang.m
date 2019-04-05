function [u1,u2,u3]=yang(sph,xi,z0,x,y,z,matrl,e_theta,coeffs,tp)
% Calculate the double force (star) and dilatation (dila) displacements U
% for a SPHEROIDAL pressure source in an elastic halfspace 
% (based on Yang et al., vol 93, JGR, 4249-4257, 1988) with arbitrary plunge (theta)
% of the long axis of the spheroid (theta = 90, prolate; theta = 0, oblate).
% Evaluate at for xi.
%
% Inputs: theta: dip angle of source
%             P: pressure change in magma chamber
%             a: semimajor axis of spheroid
%             b: semiminor axis of spheriod
%            xi: evaluate integrals at +- c
%            z0: depth of source (allowed to vary with topo)
%             x: x location of point on surface
%             y: y location of point on surface
% Output: rd: calculated range displacement
% NOTE: the x, y locations assume a source at origin
% ALSO: the units need to be in mks units so input x, y, and z0
%       in km will be changed into meters
% NOTE: In the expressions of Yang et al. the spheroid dips parallel to the y axis
%       at x=0. We will assume (initially) that it is also centered at y=0.
epsn=1e-15;

% Get spheroid information
a     = sph(1);
b     = sph(2);
c     = sph(3);
phi   = sph(4);
theta = sph(5);
Pstar = sph(6);
Pdila = sph(7);
a1    = sph(8);
b1    = sph(9);
P     = sph(10);

sinth = e_theta(1);
costh = e_theta(2);

%Poisson's ratio, Young's modulus, and the Lame coeffiecents mu and lamda
nu=matrl(3);

nu4=coeffs(2);
nu2=1-2*nu;
nu1=1-nu;
coeff=a*b^2/c^3*coeffs(1);

% Introduce new coordinates and parameters (Yang et al., 1988, page 4251):
xi2=xi*costh;
xi3=xi*sinth;
y0=0;
%z=0;
z00=tp+z0;
x1=x;
x2=y-y0;
x3=z-z00;
xbar3=z+z00;
y1=x1;
y2=x2-xi2;
y3=x3-xi3;
ybar3=xbar3+xi3;
r2=x2*sinth-x3*costh;
q2=x2*sinth+xbar3*costh;
r3=x2*costh+x3*sinth;
q3=-x2*costh+xbar3*sinth;
rbar3=r3-xi;
qbar3=q3+xi;
R1=(y1.^2+y2.^2+y3.^2).^(0.5);
R2=(y1.^2+y2.^2+ybar3.^2).^(0.5);

%C0=y0*costh+z00*sinth;  % check this!
C0=z00/sinth; 

betatop=(costh*q2+(1+sinth)*(R2+qbar3));
betabottom=costh*y1;
%atnbeta=atan2(betatop,betabottom);
%atnbeta=atan(betatop./betabottom);
nz=find(abs(betabottom)~=0);
atnbeta=pi/2*sign(betatop);
atnbeta(nz)=atan(betatop(nz)./betabottom(nz));

% Set up other parameters for dipping spheroid (Yang et al., 1988, page 4252):
% precalculate some repeatedly used natural logs:
Rr=R1+rbar3;
Rq=R2+qbar3;
Ry=R2+ybar3;
lRr=log(Rr);
lRq=log(Rq);
lRy=log(Ry);

A1star=a1./(R1.*Rr) + b1*(lRr+(r3+xi)./Rr);
Abar1star=-a1./(R2.*Rq) - b1*(lRq+(q3-xi)./Rq);
A1=xi./R1 + lRr;
Abar1=xi./R2 - lRq;
A2=R1-r3.*lRr;
Abar2=R2 - q3.*lRq;
A3=xi*rbar3./R1 + R1;
Abar3=xi*qbar3./R2 -R2;

B=xi*(xi+C0)./R2 - Abar2 - C0.*lRq;
Bstar=a1./R1 + 2*b1*A2 + coeffs(2)*(a1./R2 + 2*b1*Abar2);
F1=0;
F1star=0;
F2=0;
F2star=0;
% Skip if displacement calculated at surface (z=0)
if z ~= 0
  F1=-2*sinth*z*(xi*(xi+C0)./R2.^3 + (R2+xi+C0)./(R2.*(Rq)) + ...
     4*(1-nu)*(R2+xi)./(R2.*(Rq)));
  F1star=2*z*(costh*q2.*(a1*(2*Rq)./(R2.^3.*(Rq).^2) - ...
         b1*(R2+2*xi)./(R2.*(Rq).^2)) + ...
         sinth*(a1./R2.^3 -2*b1*(R2+xi)./(R2.*(Rq))));
  F2=-2*sinth*z*(xi*(xi+C0).*qbar3./R2.^3 + C0./R2 + (5-4*nu)*Abar1);
  F2star=2*z*(a1*ybar3./R2^3 - ...
         2*b1*(sinth*Abar1 + costh*q2.*(R2+xi)./(R2.*Rq)));
end

% calculate little f's
ff1=xi*y1./Ry + 3/(costh)^2*(y1.*lRy*sinth - ...
    y1.*lRq + 2*q2.*atnbeta) + 2*y1.*lRq - 4*xbar3.*atnbeta/costh;
ff2=xi*y2./Ry + 3/(costh)^2*(q2.*lRq*sinth - ...
    q2.*lRy + 2*y1.*atnbeta*sinth + costh*(R2-ybar3)) - ...
    2*costh*Abar2 + 2/costh*(xbar3.*lRy - q3.*lRq);
ff3=(q2.*lRq - q2.*lRy*sinth + 2*y1.*atnbeta)/costh + ...
    2*sinth*Abar2 + q3.*lRy - xi;

% Assemble into x, y, z displacements (1,2,3):
u1=coeff*(A1star + nu4*Abar1star + F1star).*y1;
u2=coeff*(sinth*(A1star.*r2+(nu4*Abar1star+F1star).*q2) + ...
       costh*(Bstar-F2star) + 2*sinth*costh*z*Abar1star);
u3=coeff*(-costh*(Abar1star.*r2+(nu4*Abar1star-F1star).*q2) + ...
       sinth*(Bstar+F2star) + 2*(costh)^2*z*Abar1star);
u1=u1+2*coeff*Pdila*((A1 + nu4*Abar1 + F1).*y1 -coeffs(3)*ff1);
u2=u2+2*coeff*Pdila*(sinth*(A1.*r2+(nu4*Abar1+F1).*q2) - ...
       coeffs(3)*ff2 + 4*nu1*costh*(A2+Abar2) + ...
       costh*(A3 - nu4*Abar3 - F2));
u3=u3+2*coeff*Pdila*(costh*(-A1.*r2 + (nu4*Abar1 + F1).*q2) + ...
       coeffs(3)*ff3 + 4*nu1*sinth*(A2+Abar2) + ...
       sinth*(A3 + nu4*Abar3 + F2 - 2*nu4*B));

