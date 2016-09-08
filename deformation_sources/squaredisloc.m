function [u]=squaredisloc(par, coord, nu)
%disloc1topo      - computes surface displacements for a dislocation source with topographic approximation 
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
%           u = matrix of displacements: Ux, Uy, Uz <length>
%
%given positive.  Keep your length units consistent! If you mix km and m you may get
%unexpected results, particularly with the strains and tilts.   
%
%21-09-2005 Falk Amelung.

%Check arguments
	%if nargin < 5 | nargin > 5 | nargout > 1 error('Usage: [u]=disloc1topo(disloc, xloc, nu, hgt, deltahgt).'); end
	%if size(par,1)~=10 error('First argument must be a 10 element dislocation geometry and slip vector.'); end
	%if size(coord,1)~=2 error('Second argument must be a 2xn matrix of station coordinates.'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       npar=zeros(1,10)';
       npar(2:10)=par(:);
       npar(1)=par(1);
       u   = disloc1(npar, coord, nu);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
