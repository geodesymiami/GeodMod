function newmat = matfill(indmat, matsize);
% matfill        - Takes in quadtree index matrix and produces full matrix
%
% function newmat = matfill(indmat, matsize);
%
% This function takes in quadtree index matrix and converts it to
% big matrix of values, to compare with the original image.
%
% INPUT:
%        indmat  - (k,m) quadtree index matrix with k pts & of level m-2
%        matsize - (1,1) size of output matrix (e.g. 1024)
%
% OUTPUT:
%        newmat  - (matsize x matsize) matrix with values from the
%                  the last column of the quadtree index matrix at 
%                  locations indicated by the m-2 first column of that
%                  same matrix


  % Get location of values, the quadtree level and the # of pts.  
  val = size(indmat,2);
  level = val - 4;
  len = size(indmat,1);  
  
  newmat = zeros(matsize)*NaN;
  
  % loop over every point
  for k=1:len
    blcksz = matsize;
    lst = 1; cst = 1;
    
    % loop over every quadtree level
    for l=1:level
      if indmat(k,l) ~= 0 
        blcksz = blcksz/2;
        switch indmat(k,l)
          case 1
	    lst = lst; cst = cst; 
          case 2	
	    lst = lst; cst = cst + blcksz;
          case 3
	    lst = lst + blcksz; cst = cst + blcksz;
          case 4
	    lst = lst + blcksz; cst = cst;
	end
      end
    end 
   
    % Fill into the new matrix
      tmp = ones(blcksz);
      [c1,c2]=find(tmp);
      d = [ones(blcksz^2,1) c1 c2]*indmat(k,val-2:val)';
      tmp(:) = d;
      newmat(lst:lst+blcksz-1,cst:cst+blcksz-1) = tmp;
    
  end
 
  











