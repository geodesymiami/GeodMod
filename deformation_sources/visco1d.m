function [u]=visco1d(par, coord, visco1dopt)
%visco1d      - computes surface displacements for a dislocation source with topographic approximation
%
%        [u]=disloc1topo(disloc, coord, nu, hgt, deltahgt)
%
%Inputs:
%         par      = dislocation parameter
%         coord    = Matrix of local station coordinates <length>, stored columnwise,
%                    i.e., east in first row, north in second row
%         nu       = Poisson's ratio
%
%Output:
%          u = matrix of displacements: Ux, Uy, Uz <length>
%
%ATTENTION INSTEAD OF VISCO1D A MOGI SOURCE IS CALLED WITH OPENING USED AS STRENGTH

nu = visco1dopt.nu ;

mogipar(1) = par(6);   % east
mogipar(2) = par(7);  % north
mogipar(3) = par(3);   % depth
mogipar(4) = par(10);  % opening is used as strength

[u] = mogi(mogipar(1:4),coord,nu);

basemap  = visco1dopt.basemap ;
coord_ll = ll2local(coord,basemap,-1);
