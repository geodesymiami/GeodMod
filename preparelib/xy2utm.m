function [md] = xy2utm(md,igram,objfuncopt,direc)
%xy2utm              - convert from xy (disloc) coordinate system into utm 
%
% INPUT: 
%       md 	             - model parameter vector (Nx1) or coordinate vector (2xN)  
%       igram            - datastruction with information on x_first, etc
%       objfunctopt      = option structure for objective function
%       dir              - direction (1: xy->utm; -1: utm->xy)          (default: 1)
%
% OUTPUT:
%       md               - model parameter or coordinate vector in UTM
%

if ~exist('direc') direc =1; end  

[x_first,y_first,x_step,y_step]=deal(igram(1).x_first,igram(1).y_first,igram(1).x_step,igram(1).y_step);
y_length = -y_step*size(igram(1).data,2) / 1000 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select which parameters are coordinates
if size(md,2)==1 
      %  Put in the fixed parameters:
      fixind=objfuncopt.fixind; fixpar=objfuncopt.fixpar; npar=zeros(length(md)+length(fixind),1); npar(fixind)=fixpar; last=1;
      for i=1:length(fixind) ind=[last:fixind(i)-1]; npar(ind)=md(ind-i+1); last=fixind(i)+1; end ; md=npar;

   [x,y]  = deal( md(6), md(7));
   [x(1),y(1)]  = deal( md(6), md(7));
   if objfuncopt.inflation   [x(2),y(2)]  = deal( md(11), md(12)); end ;
else
   [x,y]  = deal( md(1,:), md(2,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if direc==1                                       % from xy to utm
   y_length = -y_step*size(igram(1).data,2) / 1000 ; % [km]

   x_new    =  x_first    + x*1000           ;
   y_new    =  y_first    - (y_length-y)*1000;
else                                          % from utm to xy
   y_length = -y_step*size(igram(1).data,2) ;    % [m]

   x_new    =  (x-x_first)/1000              ;
   y_new    =  (y - (y_first-y_length))/1000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(md,2)==1
   [md(6),md(7)]  =deal(x_new(1),y_new(1));
   [md(11),md(12)]=deal(x_new(2),y_new(2));
else
   [md(1,:),md(2,:)]  = deal( x_new, y_new);
end
