function newindmat = check_quadtree(oldindmat, data, tolerance, fittype);
% check_quadtree - Checks if we need further quadtree partitioning
%
% function newindmat = check_quadtree(oldindmat, data, tolerance);
%  
% This script checks if we can put a constant value for
% a square instead of a quadtree it once again.  It works
% with the last two columns of the quadtree matrix, the 
% 0-1 flag control column and value column
%
% INPUT:
%            oldindmat   - (mxl) index matrix
%            data        - (nxn) data matrix
%            tolerance   - (1x1) rms tolerance
% OUPUT:
%            newindmat   - (mxl) new index matrix
  
  
  [ilin,icol] = size(oldindmat);
  
  for k=1:ilin
    if oldindmat(k,icol-3) == 0
      chunck = getchunck(oldindmat(k,:),data);
      %[c1,c2] = find(chunck==0 | chunck);   % Falk Nov 2003, this line give an error
      [c1,c2] = find(isnan(chunck)~=2);      % suggested fix from SangHo
      chunck = chunck(:);
      nn = find(isnan(chunck)==0);
      chunck_noNaN = chunck(nn);
      c1 = c1(nn); c2 = c2(nn);
      if length(chunck_noNaN) >= length(chunck)/2
	if fittype == 1 & length(chunck_noNaN) >= 3
	  [m, G, rms] = fit_bilinplane(chunck_noNaN, [c1 c2]');
	else
           medvalue = median(chunck_noNaN);
	   tmp = ones(size(chunck_noNaN)) .* medvalue;
           dif = chunck_noNaN - tmp;
           rms = sqrt( mean( dif(:).^2 ) );
	   m = [medvalue 0 0]';         
	end
        if rms <= tolerance
	  oldindmat(k,icol-3:icol) = [ 1 m' ];
	end
      elseif length(chunck_noNaN) < 1
	oldindmat(k,icol-3:icol) = [ 1 NaN NaN NaN ];
      end
    end
  end
  
  newindmat = oldindmat;





