function [U,D,S]=Mogi(m,xloc,nu,mu)
%Mogi    [U,D,S]=Mogi(m,xloc,nu,mu)
%
%Computes displacements, strains and stresses from a point (Mogi) source.
%Inputs m and xloc can be matrices; for multiple models, the deformation
%fields from each are summed.
%
%Inputs:
%    m = 4xn volume source geometry (length; length; length; length^3)
%        (x-coord, y-coord, depth(+), volume change)
% xloc = 3xs matrix of observation coordinates (length)
%   nu = Poisson's ratio
%   mu = shear modulus (if omitted, default value is unity)
%
%Outputs:
%    U = 3xs matrix of displacements (length)
%        (Ux,Uy,Uz)
%    D = 9xn matrix of displacement derivatives
%        (Dxx,Dyx,Dzx,Dxy,Dyy,Dzy,Dxz,Dyz,Dzz)
%    S = 6xn matrix of stresses
%        (Sxx,Sxy,Sxz,Syy,Syz,Szz)
%
%June 17, 1998. Peter Cervelli.
%Revised November 3, 2000.
%Fixed a bug ('*' multiplication should have been '.*'), August 22, 2001. Kaj Johnson
%Fixed error in comments that had transposed the elements of D, February 6, 2002. Peter Cervelli
%
%For information on the basis for this code see:
%
%Okada, Y. Internal deformation due to shear and tensile faults in a half-space,
%    Bull. Seismol. Soc. Am., 82, 1018-1049, 1992.

if nargin <4
     mu =1;
end

if size(xloc,1)==2
    xloc(3,:)=0;
end


%Define constant lambda

   lambda=2*mu*nu/(1-2*nu);
   U=zeros(3,size(xloc,2));
   D=zeros(9,size(xloc,2));
   S=zeros(6,size(xloc,2));

%Loop over models

     for i=1:size(m,2)

          C=m(4)/(4*pi);
          x=xloc(1,:)-m(1,i);
          y=xloc(2,:)-m(2,i);
          z=xloc(3,:);
          d1=m(3,i)-z;
          d2=m(3,i)+z;
          R12=x.^2+y.^2+d1.^2;
          R22=x.^2+y.^2+d2.^2;
          R13=R12.^1.5;
          R23=R22.^1.5;
          R15=R12.^2.5;
          R25=R22.^2.5;
          R17=R12.^3.5;
          R27=R12.^3.5;

          %Calculate displacements

               U(1,:) = U(1,:) + C*( (3 - 4*nu)*x./R13 + x./R23 + 6*d1.*x.*z./R15 );
               U(2,:) = U(2,:) + C*( (3 - 4*nu)*y./R13 + y./R23 + 6*d1.*y.*z./R15 );
               U(3,:) = U(3,:) + C*( (3 - 4*nu)*d1./R13 + d2./R23 - 2*(3*d1.^2 - R12).*z./R15);

          %Calculate their derivatives
          
          	   D(1,:) = D(1,:) + C*( (3 - 4*nu)*(-3*x.^2+R12)./R15 + (-3*x.^2+R22)./R25 + 6*z.*d1.*(-5*x.^2+R12)./R17);
               D(2,:) = D(2,:) + C*3*( (-3 + 4*nu)*(x.*y)./R15 - (x.*y)./R25 - 10*x.*y.*z.*d1./R17);
               D(3,:) = D(3,:) + C*3*( (-3 + 4*nu)*(x.*d1)./R15 - (x.*d2)./R25 - 2*x.*z.*(-5*d1.^2+R12)./R17);
               D(4,:) = D(4,:) + C*3*( (-3 + 4*nu)*(x.*y)./R15 - (x.*y)./R25 - 10*x.*y.*z.*d1./R17);
               D(5,:) = D(5,:) + C*( (3 - 4*nu)*(-3*y.^2+R12)./R15 + (-3*y.^2+R22)./R25 + 6*z.*d1.*(-5*y.^2+R12)./R17);
               D(6,:) = D(6,:) + C*3*( (-3 + 4*nu)*(y.*d1)./R15 - (y.*d2)./R25 - 2*y.*z.*(-5*d1.^2+R12)./R17);
               D(7,:) = D(7,:) + C*3*( (5 - 4*nu)*(x.*d1)./R15 - (x.*d2)./R25 - 2*x.*z.*(-5*d1.^2+R12)./R17);
               D(8,:) = D(8,:) + C*3*( (5 - 4*nu)*(y.*d1)./R15 - (y.*d2)./R25 - 2*y.*z.*(-5*d1.^2+R12)./R17);
               D(9,:) = D(9,:) + C*( (4*nu - 1)*(-3*d1.^2+R12)./R15 + (-3*d2.^2+R22)./R25 + 6*z.*d1.*(-5*d1.^2+3*R12)./R17);

          %Calculate stresses

               theta = lambda*(D(1,:)+D(5,:)+D(9,:));

               S(1,:) = S(1,:) + theta + 2*mu*D(1,:);
               S(2,:) = S(2,:) + mu*(D(2,:)+D(4,:));
               S(3,:) = S(3,:) + mu*(D(3,:)+D(7,:));
               S(4,:) = S(4,:) + theta + 2*mu*D(5,:);
               S(5,:) = S(5,:) + mu*(D(6,:)+D(8,:));
               S(6,:) = S(6,:) + theta + 2*mu*D(9,:);

   end