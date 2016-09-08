function [u,PM]=LinearPressureDike(varargin) 
%LinearPressureDike     u=LinearPressureDike(m,xy,nu,mu,nlength,nwidth)
%
%Numerically approximates the displacements for a dike opening with
%constant pressure. 'm' specifies the dike geometry and pressure,
%'xy' is a 2xn matrix of observation coordinates, 'nu' is Poisson's
%ratio, 'mu' is the shear modulus, and 'n' defines the precision
%of the numerical approximation.  The last four arguments can be
%omitted; they default to 0.25, 3e10 Pa, 4,4 respectively.
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
%  m(9) = PressureGradient  (force / length^2)
%
%Warning!  Keep length units consistent among the dike parameters and
%the observation coordinates.
%  FA April 2002 based on code from Peter Cervelli
%  FA June  2002, modified such that it works with Dike dimensions given in km
%           Unit problem: This routine calls disloc3d which needs dislocation dimension
%           and slip vector in same units (all in km orall in meter).
%           It also calls disloc which uses slipvector (opening) in meter but dimensions in km.
%           Cervelli uses originall DM_Disloc which probably (?) uses disloc dimension and slipvector in same unit
%

m=varargin{1};
xy=varargin{2};
varargin{7}=[];

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
   nlength=4;
else
   nlength=varargin{5};
end

if isempty(varargin{6})
   nwidth=4;
else
   nwidth=varargin{6};
end
%s=sprintf( 'm: %4.1f %4.1f %4.1f %6.1f %6.1f %4.1f %4.1f %6.2e ',m(1:8)); disp(s)
%Define dike grid

   nhat=unitvector(m);
   PM=patchfault(m,nlength,nwidth);
   PM(8:9,:)=0;
   PM(10,:)=1;
   x=midpoint(PM);
   x(3,:)=-x(3,:);

%Make kernel for constant pressure

   tmpPM=PM ; tmpPM([1 2 3 6 7],:)=tmpPM([1 2 3 6 7],:)*1000; x=x*1000;
   for i=nlength*nwidth:-1:1
      [U,D,sigma,flag]=disloc3d(tmpPM(:,i),x,mu,0.25);
%      [U,D,sigma,flag]=disloc3d(PM(:,i),x,mu,0.25);
      t=[nhat'*sigma([1 2 3],:)
         nhat'*sigma([2 4 5],:)
         nhat'*sigma([3 5 6],:)];
      S(:,i)=(nhat'*t)';
   end

%Solve for opening distribution

   del_width=m(2)/nwidth;
   %press_grad=[0:nwidth-1] * del_width*m(9);
   press=m(8)+[0:nwidth-1] * del_width*m(9);
   press=repmat(press',nlength,1);
   op=S\press;
%   disp([min(op);max(op)]');
   PM(10,:)=op';

   u=DM_Disloc(PM,xy,nu);
