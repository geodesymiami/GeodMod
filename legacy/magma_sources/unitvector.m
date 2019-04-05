function u=unitvector(m)
%UNITVECTOR     u=unitvector(m)
%
%Returns a unit vector perpendicular to the plane of
%dislocation m.  Models should be stored columnwise
%for vectorized output.

%Parse inputs 

   if nargin ~=1
      help unitvector
   end

   [i,j]=size(m);

   if j>i & min([i j])==1
      m=m';
   end
       
%Calculate unit vector

   theta=(90-m(5,:))*pi/180;
   delta=m(4,:)*pi/180;
   sd=sin(delta);
   cd=cos(delta);
   st=sin(theta);
   ct=cos(theta);
   cd(abs(cd)<1e-15)=0;
   sd(abs(sd)<1e-15)=0;
   ct(abs(ct)<1e-15)=0;
   st(abs(st)<1e-15)=0;
   u(1,:)=st.*sd;
   u(2,:)=-ct.*sd;
   u(3,:)=cd;