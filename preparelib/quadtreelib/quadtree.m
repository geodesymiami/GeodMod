function [indexmatrix] = quadtree(oldind)
% quadtree       - Add a new quadtree partitioning level 
%
% function indexmatrix = quadtree(oldind)  
%
% Input:    
%           oldind      - (m x l) old index matrix of level l-2
% 
% Output:
%           indexmatrix - (t x l+1) new index matrix of level l-1
%                         where t>m.
		
  indexmatrix = [];
  [lin col] = size(oldind);
  nlin = 1;
  
  % loop over every old quadtree partition
  for k = 1:lin
    if oldind(k,col-3) == 1  % If deeper part. isn't needed, we add a 0
      tmp1 = [oldind(k,1:col-4), 0];
      tmp2 = oldind(k,col-3:col);
      
      indexmatrix = [indexmatrix ; tmp1 tmp2];
      nlin = nlin+1;
    else % Deeper partition needed, we add three new lines to the matrix
      tmp1 = [ones(4,1)*oldind(k,1:col-4), [1 2 3 4]'];
      tmp2 = [zeros(4,1) ones(4,1)*oldind(k,col-2:col)];
      
      indexmatrix = [indexmatrix ; tmp1 tmp2];
      nlin = nlin+4;
    end
  end
  

