	function u = abs_disp(nu, dis_geom, xy)

%  modified in April 2002 by FA: eliminating for loop.  Only checked for opening dislocation

	[nsta,m] = size(xy);
 
	u = zeros(3,nsta);
	u = disloc1(dis_geom', xy', nu); 

%	for i= 1:nsta
%		u(:,i) = disloc(dis_geom', [xy(i,1); xy(i,2)], nu); 
%	end


%	u = u';
