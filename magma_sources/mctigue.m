function U = mctigue(m,obs,mu,nu)

%MCTIGUE Calculate surface deformation due to a pressurized finite sphere.
%
%   U = mctigue(m,obs,nu) calculate surface displacement field at 
%   observation points due to a pressurized finite sphere.  The result is
%   similar to Mogi source, but this function include higher order terms.
%
%INPUT:
%   m(1)    east
%   m(2)    north
%   m(3)    depth(+)
%   m(4)    radius
%   m(5)    excess pressure, same unit as mu
%
%   obs     2 x number of observation points
%   mu      shear modulus
%   nu      Poisson's ratio
%
%OUTPUT:
%   U       surface displacement, 3 x size(obs,2)
%
%
% Sang-Ho Yun, September 29, 2005
%
%
% For information on the basis for this code see:
%
% MCTIGUE, D. F., 1987. Elastic stress and deformation near a finite 
% spherical magma body: resolution of the point source paradox. JGR 
% 92(B12), p.12931-40.
%
% Williams, C. A., Wadge, G., 1998. The effects of topography on magma
% chamber deformation models: Application to Mt. Etna and radar 
% interferometry. GRL, 25(10), p.1549-52. 

if (nargin == 0)
    help mctigue;
    return;
end

if size(obs,1)~=2
    error('The dimension of obs should be 2 x n');
end

east = m(1);
north = m(2);
d = m(3);
a = m(4);
p = m(5);

x = obs(1,:) - east;
y = obs(2,:) - north;
[th,r] = cart2pol(x,y);
rho = r/d;
aod = a/d;

aod3 = aod^3;
aod6 = aod^6;
denom1 = rho.^2 + 1;
denom2 = 7-5*nu;

common1 = p*(1-nu)./(mu*denom1.^1.5);
common2 = aod3 - 0.5*aod6*(1+nu)/denom2 + 3.75*aod6*(2-nu)/denom2./denom1;
Uz = d*common1.*common2;
Ur = r.*common1.*common2;

[Ux,Uy] = pol2cart(th,Ur);

U = [Ux; Uy; Uz];