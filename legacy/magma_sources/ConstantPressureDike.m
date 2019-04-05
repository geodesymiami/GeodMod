function [u,PM]=ConstantPressureDike(varargin) 
%ConstantPressureDike     u=ConstantPressureDike(m,xy,nu,mu,n)
%
%Numerically approximates the displacements for a dike opening with
%constant pressure. 'm' specifies the dike geometry and pressure,
%'xy' is a 2xn matrix of observation coordinates, 'nu' is Poisson's
%ratio, 'mu' is the shear modulus, and 'n' defines the precision
%of the numerical approximation.  The last three arguments can be
%omitted; they default to 0.25, 3e10 Pa, and 4 respectively.
%
%The dike is defined as follows:
%
%  m(1) = Length    (length)
%  m(2) = Width     (length)
%  m(3) = Depth     (length)
%  m(4) = Dip       (degrees)
%  m(5) = Strike    (degrees)
%  m(6) = Easting   (length)
%  m(7) = Northing  (length)
%  m(8) = Pressure  (force / length^2)
%
%Warning!  Keep length units consistent among the dike parameters and
%the observation coordinates.

m=varargin{1};
xy=varargin{2};
varargin{6}=[];

if isempty(varargin{3})
   nu=0.25;
else
   nu=varargin{3};
end

if isempty(varargin{4})
   mu=3e10;
else
   mu=varargin{4};
end

if isempty(varargin{5})
   n=3;
else
   n=varargin{5};
end

%s=sprintf( 'm: %4.1f %4.1f %4.1f %6.1f %6.1f %4.1f %4.1f %6.2e ',m(1:8)); disp(s)
%Define dike grid

   nhat=unitvector(m);
   j=2^n;
   PM=patchfault(m,j,j);
   PM(8:9,:)=0;
   PM(10,:)=1;
   x=midpoint(PM);
   x(3,:)=-x(3,:);

%Make kernel for constant pressure

   for i=j^2:-1:1
      [U,D,sigma,flag]=disloc3d(PM(:,i),x,mu,0.25);
      t=[nhat'*sigma([1 2 3],:)
         nhat'*sigma([2 4 5],:)
         nhat'*sigma([3 5 6],:)];
      S(:,i)=(nhat'*t)';
   end

%Solve for slip distribution

   s=S\repmat(m(8),j^2,1);
   max(s);
   PM(10,:)=s';

   u=DM_Disloc(PM,xy,nu);
