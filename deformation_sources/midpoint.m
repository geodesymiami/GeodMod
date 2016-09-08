function [mp]=midpoint(m)
%MIDPOINT    mp=MIDPOINT(m)
%
%Returns coordinates at the mid-point of dislocation m.
%Models should be stored columnwise for vectorized output.

%Parse inputs 

   if nargin ~=1
      help midpoint
   end

   [i,j]=size(m);

   if j>i & min([i j])==1
      m=m';
   end
       
%Calculate midpoints

   theta=(90-m(5,:))*pi/180;
   delta=m(4,:)*pi/180;
   W=m(2,:)/2;
   mp(1,:)=m(6,:)-W.*cos(delta).*sin(theta);
   mp(2,:)=m(7,:)+W.*cos(delta).*cos(theta);
   mp(3,:)=m(3,:)-W.*sin(delta);



