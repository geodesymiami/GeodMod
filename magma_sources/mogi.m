function [u,e,t]=Mogi(volgeom, xloc, nu)
%Mogi    [u,e,t]=Mogi(volgeom, xloc, nu)
%
%Computes surface displacements, strains, and tilts due to a Mogi source.
%
%Inputs:
%     volgeom = Mogi source geometry: East, North, Depth, Volume change
%               <length, length, length, volume>
%        xloc = Matrix of local station coordinates <length>, stored columnwise,
%               i.e., east in first row, north in second row
%          nu = Poisson's ratio
%
%Output:
%           u = matrix of displacements: Ux, Uy, Uz <length>
%           e = matrix of strains: Exx, Exy, Eyy
%           t = matrix of tilts: dUz/dx, dUz/dy
%
%Notes: The term 'depth' denotes an unsigned length and should therefore always be
%given positive.  Keep your length units consistent! If you mix km and m you may get
%unexpected results, particularly with the strains and tilts.   
%
%10-18-00 Fixed a bug in the tilts. PFC
%05-17-98 Peter Cervelli.

%Check arguments

	if nargin < 3 | nargin > 5 | nargout > 3
		error('Usage: [u,e,t]=Mogi(volgeom, xloc, nu).');
	end

	[v]=size(volgeom);

	if min(v)~=1 | max(v)~=4
		error('First argument must be a 4 element source geometry vector.');
	end

	[x]=size(xloc);
	
	if x(1)~=2
		error('Second argument must be a 2xn matrix of station coordinates.');
	end

	[P]=length(nu);

	if P~=1
		error('Third argument must be a scalar (Poisson''s ratio).');
	end  	

%Compute displacements

	E=volgeom(1)-xloc(1,:);
	N=volgeom(2)-xloc(2,:);
	E2=E.^2;
	N2=N.^2;
	d2=volgeom(3)^2;
	C=(nu-1)*volgeom(4)/pi;
	R=sqrt(d2+E2+N2);
	R3=C*R.^-3;
	u=[E.*R3;N.*R3;-volgeom(3)*R3];

%Compute strains (if necessary)

	if nargout > 1
 		R5=C*R.^-5;
     	e(1,:)=R5.*(2*E2-N2-d2);
     	e(2,:)=3*R5.*(E.*N);
     	e(3,:)=R5.*(2*N2-E2-d2);
	end

%Compute tilts (if necessary)

	if nargout > 2
   	t(1,:)=-3*R5.*E*volgeom(3);
   	t(2,:)=-3*R5.*N*volgeom(3);
	end
