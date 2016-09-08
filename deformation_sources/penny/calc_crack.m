% Matlab routine for calculating surface displacements due to a
% uniformly pressurized horizontal penny-shaped crack
% for details, see the 2001 GJI paper.
% v. 1.0 Yuri Fialko 7/11/2000

% The following parameters need to be user-supplied:   
nis=2;       % number of sub-intervals on [0,1] on which integration is done
             % using a 16-point Gauss quadrature (i.e., total of nis*16 points)
eps=1e-5;    % solution accuracy for Fredholm integral equations (stop 
             % iterations when relative change is less than eps)
h=1;         % dimensionless crack depth (Depth/Radius ratio)
r=[0:0.1:3]; % observation points at the surface 

% Solve a coupled system of Fredholm eqs. of 2nd kind for basis functions fi,psi
% t and Wt are nodes and weights of the num. integration quadrature
 [fi,psi,t,Wt]=fredholm(h,nis,eps);

 [Uv,Ur]=intgr(r,fi,psi,h,Wt,t);  %calculate vertical and radial displ.

