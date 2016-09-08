function [dataset]=qtsjonni2qtgeodmod(inData,dataset, basemap, location)
%qtsjonni2qtgeodmod  - converts Sjonni's quadtree coordinates into geodmod's coordinate system
%
%  FA Jan 2011    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin==3) location='Haiti'; end

switch location

    case {'Haiti'}
      dataset.datavec   = inData.sv';
      dataset.radarlook = inData.los;
      dataset.cov       = inData.cov;
      dataset.Ndata     = length(dataset.datavec);

      dataset.coord     = [round(inData.cnt(1,:)/basemap.x_posting); round(-inData.cnt(2,:)/basemap.y_posting)] ;
      dataset.cx        = round(inData.cx/basemap.x_posting);
      dataset.cy        = round(-inData.cy/basemap.y_posting);
      %dataset.cy        = round((inData.cy+size(basemap.data,1)*basemap.y_posting)/basemap.y_posting);
    end
