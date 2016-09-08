function [Lap, Lap_inv] = modelwt2(nve, nhe, delx, dely, surf)
%%
%%	 [Lap, Lap_inv] = modelwt(nve, nhe, delx, dely, surf);
%%
%%  Model weighting matrix to compute Laplacian and its inverse
%%  for smoothing in slip inversions.
%%
%%  Input:
%%		nve 	= number of vertical elements
%%		nhe 	= number of horizontal elements
%%		delx 	= length of elements in along strike dimension
%%		dely 	= length of elements in dip dimension
%%		surf 	= 1 if fault breaks free surface ~= 1 otherwise
%%  Output:
%%		Lap 	= finite difference Laplacian in two dimensions
%%		Lap_inv = inverse of Lap

	ngrid = nhe * nve;
	Lap = zeros(ngrid);  xpartd = zeros(ngrid); ypartd = zeros(ngrid);
	
%       x-derivative for central part of grid
%       (exclude left & right edges)
        for i  = nve+1:(nhe*nve-nve)
		xpartd(i,i-nve) = + 1.0;
		xpartd(i,i)     = - 2.0;
		xpartd(i,i+nve) = + 1.0;
 	end

%       y-derivative for central part of grid
%       (exclude top & bottom edges)
	temp = zeros(nve);
	for i = 2:nve-1
		temp(i,i-1:i+1) = [1 -2 1];
	end
	for j = 1:nhe
		k = (j-1)*nve; 
		ypartd( k+1:k+nve, k+1:k+nve ) = temp;
	end

%% 	need to do top edge y-derivative
%%	 for slip breaking the surface
   for i  = 1:nhe
		k = i*nve;
		if surf == 1 
			ypartd( k, k-1:k) = [1 -1];
		else
			ypartd( k, k-1:k) = [1 -2];
		end
 	end


%%	bottom edge y-derivative
        for i  = 1:nhe
		k = (i-1)*nve+1; 
		ypartd( k, k:k+1) = [-2 1];
 	end


%%	left edge x-derivative
        for i  = 1:nve
		xpartd(i,i) 	    = - 2.0;
		xpartd(i,i+nve)     = + 1.0;
 	end
	
%%	right edge x-derivative
        for i  = (nhe-1)*nve+1:nhe*nve
		xpartd(i,i) 	    = - 2.0;
		xpartd(i,i-nve)     = + 1.0;
 	end
    
  % Check if elements have various lengths
  if length(delx) > 1
     %delxx = ones(nve,1)*delx'.^2;
     xpartdd = xpartd./( delx.^2 * ones(1,nve*nhe) );
%     xpartdd = xpartd./( ones(nve*nhe,1) * delx'.^2 );
  else
     xpartdd = xpartd/delx^2;
  end
  
  % Check if elements have various lengths
  if length(dely) > 1
     %xdelyy = ones(nhe,1)*dely'.^2;
     ypartdd = ypartd./( dely.^2 * ones(1,nve*nhe) );
  else
     ypartdd = ypartd/dely^2;
  end

   Lap = xpartdd + ypartdd;
   %Lap = xpartd/.(ones(delx^2 + ypartd/dely^2;
   Lap_inv = inv(Lap);
