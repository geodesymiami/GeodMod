function [xy_out] = lola2xy(xy,data,direc)
%lola2xy           - convert from lola coordinates into (local, disloc) 
%
% INPUT: 
%       xy 	             - coordinate vector (2xN)  
%       igram            - datastructure with information on x_first, etc
%       dir              - direction (1: lola->xy; -1: xy->lola)          (default: 1)
%
% OUTPUT:
%       xy_out           - model parameter or coordinate vector in UTM
%
% FA, May 1 2007                                                                                                                             
% TODO: Another way would be to first convert LL to UTM and generate local coordinates by 
% referring (subtracting) to  the lower-right corner

if ~exist('direc') direc =1; end  

   [x_first,y_first,x_step,y_step]= deal(data(1).x_first,data(1).y_first,data(1).x_step,data(1).y_step);
   [x_posting,y_posting]          = deal(data(1).x_posting,data(1).y_posting);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Following the simplified lola2xy (old: ll2local) conversion in Dem2Basemap.m, the conversion is done using the
%  distance [km] in x (long) direction and y (lat direction at the center of the area give in 
%  the data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   x_last     = x_first + x_step*size(data(1).data,2);
   y_last     = y_first + y_step*size(data(1).data,1);
   x_unitdist = 1/x_step*x_posting ;    
   y_unitdist =-1/y_step*y_posting ;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if direc==1                                       % from lola to xy 
   x_new      = (xy(1,:)-x_first)*x_unitdist ;
   y_new      = (xy(2,:)-y_last)*y_unitdist  ;
   xy_out     = [x_new;y_new];
elseif direc==-1                                       % from xy to lola %%Not yet implemented
   x_new      = x_first+xy(1,:)/x_unitdist ;
   y_new      = y_last +xy(2,:)/y_unitdist ;  
   xy_out     = [x_new;y_new];
end
