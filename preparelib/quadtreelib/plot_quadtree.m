function [coord, coordl,cx,cy] = plot_quadtree(indexmatrix,data);
% plot_quadtree  - Gives coordinates of the quadtree squares 
%
% function [coord,coordl] = plot_quadtree(indexmatrix, data);
%
% INPUT
%         indexmatrix - (1xl) quadtree index
%         data        - (nxn) data matrix, n = 2^p
% 
% OUPUT
%         coord       - (lx2) quadtree obs. locations
%         coordl      - (6l x 2 ) quadtree square outlines
%                  
  
  coord = []; coordl = []; cx = []; cy = [];
  [lin,col] = size(data);
  
  % number of points and level of quadtree partitioning
  len = size(indexmatrix,1);
  level = size(indexmatrix,2) - 4;
  
  % Loop over each quadtree point
  for k=1:len
    blcksz = lin;
    lst = 1; cst = 1;
    
    % Loop over each quadtree level
    for l = 1:level
      if indexmatrix(k,l) ~= 0    % if not zero, we go deeper
        blcksz = blcksz/2;
        switch indexmatrix(k,l)
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
   coord = [ coord ; lst - 1 + blcksz/2, cst - 1 + blcksz/2 ];
   coordl = [coordl ; NaN, NaN; lst-1, cst-1 ; lst-1, cst-1+blcksz; 
             lst-1+blcksz, cst-1+blcksz; lst-1+blcksz, cst-1; lst-1, cst-1];
   cx = [cx, [cst-1; cst-1+blcksz; cst-1+blcksz; cst-1]];
   cy = [cy, [lst-1; lst-1; lst-1+blcksz; lst-1+blcksz]];
 end
 
 
 % Plot quadtree partitioning figures:
 figure
 plot(coord(:,2),-coord(:,1),'.');
 title('Obs. locations after quadtree partioning')
 axis image
 axis([0 col -lin 0])
 
 
 figure
 plot(coordl(:,2),-coordl(:,1));
 title('Quadtree partitioning')
 axis image
 axis([0 col -lin 0])




