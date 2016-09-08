function [u]=yang_sourcetopo(par, coord, nu, hgt, toporanges)
%disloc1topo      - computes surface displacements for a dislocation source with topographic approximation 
%
%        [u]=disloc1topo(disloc, coord, nu, hgt, deltahgt)
%
%Inputs:
%         par      = dislocation parameter
%         coord    = Matrix of local station coordinates <length>, stored columnwise,
%                    i.e., east in first row, north in second row
%         nu       = Poisson's ratio
%         hgt      = station height                         [km]
%         toporanges= topo depth ranges  
%
%Output:
%           u = matrix of displacements: Ux, Uy, Uz <length>
%
%          Source depth counts from the lowest value in hgt
%
%
%given positive.  Keep your length units consistent! If you mix km and m you may get
%unexpected results, particularly with the strains and tilts.   
%
%21-09-2005 Falk Amelung.

%Check arguments
	if nargin < 5 | nargin > 5 | nargout > 1 error('Usage: [u]=yang_sourcetopo(disloc, xloc, nu, hgt, deltahgt).'); end
	if size(par,1)~=8 error('First argument must be a 10 element dislocation geometry and slip vector.'); end
	if size(coord,1)~=2 error('Second argument must be a 2xn matrix of station coordinates.'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    deltahgt=(max(hgt)+0.001-min(hgt))/toporanges;
    meanhgt=mean(hgt);
    ioff=(meanhgt-min(hgt))/deltahgt;
    u=zeros(3,length(coord));

    for i=-ioff:-ioff+toporanges-1
       delhgtmin  = meanhgt+i*deltahgt;
       delhgtmax  = delhgtmin+deltahgt;
       ind     = find(hgt>=delhgtmin & hgt <delhgtmax);
       
       ndep = par(3)+mean(hgt(ind))-meanhgt;
       npar =[par(1:2)', ndep, par(4:8)']';
       ncoord     = coord(:,ind) ;
       u(:,ind)   = yang_source(npar, coord(:,ind), nu);
      %[i ndep mean(hgt(ind)) delhgtmin delhgtmax length(ind)]
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
