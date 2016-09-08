function [u] = multidisloc(multidisloc_par, coord, multidislocopt)
%multidisloc    - multiple dislocations 
%
%Computes surface displacements due to 2 dislocation sources
%
%Inputs:
%Output:
%           u = matrix of displacements: Ux, Uy, Uz <length>
%
% 10/2008  Falk Amelung

 x_unit       = 'km';
 nu           = 0.25;
 N_disloc     = multidislocopt.N_disloc;
 u_third      = 0;
 
 [disloc_par] = multidislocpar2dislocpar(multidisloc_par,multidislocopt,x_unit);

 u_first  = disloc1(disloc_par(1:10), coord,nu);
 u_second = disloc1(disloc_par(11:20),coord,nu);
 
 if (N_disloc==3) u_third = disloc1(disloc_par(21:30),coord,nu); end

 u        =  u_first + u_second + u_third;
