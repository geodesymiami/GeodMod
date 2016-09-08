function [newindmat,sqval,cx,cy,cntp,nlin] = quadtree_part(data,tolerance,fittype,startlevel,maxdim);
% quadtree_part  - Quadtree partitioning of data with given tolerance
%
% function [newindmat,cc,cx,cy] = quadtree_part(data,tolerance);
%
% INPUT:
%          data      - (MxN) data matrix (can have NaN's)
%          tolerance - (1x1) tolerance value for each square
%          fittype   - (1x1) type of fit (constant=0; bilinear=1)
%          startlevel- (1x1) starting quadtree level (default =1)
%
% OUTPUT:
%          newindmat - (txs) quadtree index matrix of level s-2
%                      and t squares, the last column gives the 
%                      values for each square.
%          sqval     - (1xt) row vector of quadtree squares values
%          cx        - (4xt) x-coordinate of square corners
%          cy        - (4xt) y-coordinate of square corners
%          cntp      - (tx2) center point coordinates of squares
%          nlin      - (1x1) # number of lines in working image (e.g. 512)
%
% This function need function 'check', 'quadtree', 'plot_quadtree', and
% 'get_chunck'.
%
% SJ, 21 apr, 2000
% FA Sep 2006, changed disp by logmessage
%              commented out display
%          

if (nargin==3 | nargin==4 | nargin==5); else,help quadtree_part;return;end;
if nargin==3,startlevel=1;maxdim=13;end
if nargin==4,maxdim=13;end

%maxdim = 13;  % don't expect more than 2^13 ~ 8000 pixels a side

[lin,col]=size(data);
logmessage(' ')
logmessage(sprintf(' Size of data,      lines and columns: %d %d',lin,col))

% Here we add NaN's to change data size to (2^k x 2^k)
dim=1; condition = max([lin col]);
while condition > 2^dim 
    dim = dim+1;
end

nlin = 2^dim;
ncol = nlin;

dataexp = zeros(nlin,ncol)*NaN;
dataexp(1:lin,1:col) = data;

logmessage(sprintf(' Data working size, lines and columns: %d %d',nlin,ncol))

% This is the starting quadtree index matrix (values 
% 10 are arbitrary).  The first number tells you the
% quadrant (clockwise) and the second last value is
% a test index (zero = further partioning needed), and
% the last three columns are the bilinear values a, b, 
% and c in  a + bx + cy

indmat = [1 0 10 0 0; 2 0 10 0 0; 3 0 10 0 0; 4 0 10 0 0]; 

% Here we add to the index matrix in case don't want
% to start at the top level.
if startlevel>1
  for k = 2:startlevel
    nindmat = quadtree(indmat);
    indmat = nindmat;
  end 
end

% Here we loop over each k in 2^k.
for k = startlevel:maxdim
  % First we check if further partioning is needed and 
  % pick out check column (che)
  newindmat = check_quadtree(indmat, dataexp, tolerance, fittype);
  che = newindmat(:,size(indmat,2)-3);  
  % If any zeros in check-column then perform
  if prod(che) == 0
    % Increase by one partioning level
    indmat = quadtree(newindmat);    
  else
    k = maxdim;
  end
end


% Here we plot the quadtree points and squares:
[cntp,co2,cx,cy] = plot_quadtree(newindmat,dataexp);

% Here we plot everything with patch:
sqval = newindmat(:,size(newindmat,2)-2)';
%falk figure ; patch(cx,-cy,sqval); axis image





  

  

    


