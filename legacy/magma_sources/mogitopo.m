function [u]=mogitopo(par, coord, nu, hgt, toporanges)
%mogitopo      - computes surface displacements for a mogi source including topographic approximation
%
%Mogi    [u]=disloc1topo(disloc, c, nu, hgt, toporanges)
%
%Inputs:
%         par      = Mogi source geometry: East, North, Depth, Volume change
%                    <length, length, length, volume>
%         coord    = Matrix of local station coordinates <length>, stored columnwise,
%                    i.e., east in first row, north in second row
%         nu       = Poisson's ratio
%         hgt      = station height                         [km]
%         toporanges= height range to summarize data points  
%
%Output:
%          u = matrix of displacements: Ux, Uy, Uz <length>
%
%          Source depth counts from the lowest value in hgt
%
%given positive.  Keep your length units consistent! If you mix km and m you may get
%unexpected results, particularly with the strains and tilts.
%
%21-09-2005 Falk Amelung.

%Check arguments
	if nargin < 5 | nargin > 5 | nargout > 1 error('Usage: [u]=mogitopo(disloc, xloc, nu, hgt, deltahgt).'); end
	if size(par,1)~=4 error('First argument must be a 4 element Mogi source geometry and strength.'); end
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
       npar =[par(1:2)', ndep, par(4)']';
       ncoord     = coord(:,ind) ;
       u(:,ind)   = mogi(npar, coord(:,ind), nu);
       %[i ndep mean(hgt(ind)) delhgtmin delhgtmax length(ind)]
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
