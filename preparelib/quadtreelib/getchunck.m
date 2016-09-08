function chunck = getchunck(index, data);
% getchunck      - Gets a piece of data according to the quadtree index
%
% function chunck = getchunck(index, data);
%
% INPUT
%         index = (1xl) quadtree index
%         data  = (nxn) data matrix
% 
% OUPUT
%         chunck = (kxk) chunk as indicated by 
%                  the quadtree index
%
  
  len = length(index) - 4;
  [lin,col] = size(data);
  blcksz = lin;
  lst = 1; cst = 1;
  
  for k=1:len
    blcksz = blcksz/2;
    switch index(k)
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
  
  chunck = data(lst:lst+blcksz-1,cst:cst+blcksz-1);
  











