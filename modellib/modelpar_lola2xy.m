function [par_out] = modelpar_lola2xy(par,basemap,modelopt,direc)
%modelpar_lola2xy   - converts model parameter from xy (local, disloc) coordinate system into lola and vice-versa
%
% INPUT: 
%       par 	  - model parameter vector (Nx1) or coordinate vector (2xN)  
%       basemap   - structure with information on x_first, etc
%       modelopt  - structure with modeloptions
%       direc     - direction (1: lola->xy; -1: xy->lola)          (default: 1)
%
% OUTPUT:
%       par_out   - model parameter with converted coordinates
%
% FA, June  2007                                                                                                                             

ind          = find_coord_index(modelopt) ;

par_out      = par;
tmp_in       = zeros(2,length(ind)/2);
tmp_in(:)    = par(ind);
tmp_lola     = lola2xy(tmp_in,basemap,direc);
par_out(ind) = tmp_lola;
