function [varargout]  =  resampledata(resamplefac,varargin)
%resampledata              - resamples data on smaller or lower grid
%
% INPUT: 
%       resamplefac 	 - factor for new grid (0.2 is 5 times smaller and 2 is twicw bigger) 
%       data             - can be a regular 2D dataset, an rgb image ((:,:,3) or jpg image
%
% OUTPUT:
%       ndata1,ndata2,...  - resampled data
%
% FA, 7/2008                                                                                                                             

  data = varargin{1};
  xi           = linspace(1,size(data,2),size(data,2)*resamplefac);
  yi           = linspace(1,size(data,1),size(data,1)*resamplefac);
     %x        = linspace(1,size(data,2),size(data,2));     % not necessary becaue interp2 assumes equally spaced 
     %y        = linspace(1,size(data,1),size(data,1));     % meshgrid (plaid) format
     %[xx,yy]  = meshgrid(x,y);
     %[xxi,yyi]= meshgrid(xi,yi);

  for i=1:length(varargin)  
      data            = varargin{i};
      
      if size(data,1)>1                                   % for regular 2-D data array (or 2-D-by-3 rgbim)
         ndata           = zeros(length(xi),length(yi),size(data,3));
         for j=1:size(data,3)                             % size(data,3) is 3 for rgb image
             ndata(:,:,j)= interp2(data(:,:,j),xi',yi);       
         end
      else
             ndata       = interp1(data,xi');             % for 1-D vecor (x,y-axis)
      end

      varargout(i)    = {ndata};
  end

