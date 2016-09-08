function [m, G, rootms] = fit_bilinplane(data, coord)
% fit_bilinplane - Fits bi-linear plane to a given image
%
% function [m, G, rootms] = fit_bilinplane(data, coord)
%
% Input:
%		data   - (nx1) vector of data
%		coord  - (2xn) matrix of coordinate locations
%
% Output:		
%		m      - (3x1) vector of model parameters
%			 [a b c]' for plane z(x,y) = a + bx + cy 
%		G      - (tx3) constuction matrix with t <= n lines
%                        corresponding to no-NaN data
%               rootms - (1,1) root mean square error between plane and data

% Clean data and coordinates of NaNs
NoNanInd=find(isnan(data)==0);
d=data(NoNanInd);
coord=coord(:,NoNanInd);

N = length(d);

% If 3 or more data points left after NaN screening, then
if N >= 3
   % Make matrix G
   G = [ones(N,1) coord(1,:)' coord(2,:)' ];
   gtginv = inv(G'*G);

   m = gtginv*(G'*d);

   % Calculate the rms
   rootms = sqrt( mean ( (d-G*m).^2 ) );
   
   % Calculate the median value
   %valint=median(d);
   %corint=mean(coord');
else    
   rootms = 0; G=0; m=[0 0 0]';
end









