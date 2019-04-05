function [u] = penny_crack(param,coord)
%   penny_crack   - displacement due to a penny-shaped crack            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [u] = penny_crack(param,coord)
%
% Surface Displacements Due to Penny-Shaped Crack 
% FA, Mar 2002. Using code of Yuri Fialko (Fialko et al., GJI 2001)
% FA, May 2005. Changes do that the last model parameter is the linear one (strength)
%
%  Input: 
%        param    4*1 parameters for mogi source
%                     param(1) :  x-coordinate of source (km)
%                     param(2) :  y-coordinate of  source (km)
%                     param(3) :  depth source (km)
%                     param(4) :  radius of source
%                     param(5) :  strength of source
%        coord    2*N array  with x,y coordinates of datavec 
%                     where the model displacement will be computed
%  Output:
%        u        3*N array with displacement coordinates   [ux uy uz]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nis=2; eps=1e-5;  r=[0.001:0.1:3];

x=param(1) ;
y=param(2) ;
h=param(3) ;  % depth of crack
r=param(4) ;  % radius of crack
f=param(5) ;  % strength
xvec=coord(1,:);
yvec=coord(2,:);

d=sqrt( (xvec-x).^2 + (yvec-y).^2 );
theta=atan2(yvec-y,xvec-x);

h=h/r; d=d/r;
[fi,psi,t,Wt]=fredholm(h,nis,eps); [Uv,Ur]=intgr(d,fi,psi,h,Wt,t);

ux = Ur.*cos(theta);
uy = Ur.*sin(theta);

u=[ux ; uy; -Uv].*f;
