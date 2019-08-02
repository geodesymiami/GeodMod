function [enu,coord,u]=ForwardModel_forBasemap(dataset,modelopt,basemap,GridRowsCols,x_unit)
%function [enu,coord,u]=ForwardModel_forBasemap(modelopt,basemap,GridRowsCols,x_unit,dataset)  %%Added dataset Anieri 4/27/15
%ForwardModel_forBasemap -  plots data, model and residual in xy (local) coordinate system (using PlotModel)
%
%usage:  [enu,coord,u]=ForwardModel_forBasemap(modelopt,basemap,GridRowsCols)
%
%        modelopt     model structure
%        inverseopt   is only used for some calculations in PlotModel. Should be eliminated.
%        dataset      data structure
%
% Falk Amelung, May 2007
  logmessage(sprintf('[]=%s(%s)',mfilename,inputname(1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
  CoordSystem='Disloc';
  if nargin<4 |  isempty(GridRowsCols)
     GridRowsCols=[40 40] ;      % default
     GridRowsCols=fliplr(round(size(basemap.data)/8)) ;
  end       
  GridRows = GridRowsCols(1);
  GridCols = GridRowsCols(2);

  x_length = size(basemap.data,2)*basemap.x_posting;
  y_length = size(basemap.data,1)*basemap.y_posting;
  x        = linspace(0,x_length,GridCols);
  y        = linspace(0,y_length,GridRows);
  [xx,yy]  = meshgrid(x,y);
  xi       = linspace(0,x_length,size(basemap.data,2));
  yi       = linspace(0,y_length,size(basemap.data,1));
  [xxi,yyi]= meshgrid(xi,yi);
  coord    = [xx(:) yy(:)]';

  %hgt = dataset.hgt; %%%Anieri 4/27/15
  hgt = basemap.data; %%% FA 7/2019
  %%%
  %[~,~,u] = ForwardModel(modelopt.par.xy,coord,ones(3,length(coord)),modelopt);
  [~,~,u] = ForwardModel(modelopt.par.xy,coord,ones(3,length(coord)),modelopt,hgt); %%%Added hgt %Anieri 4/27/15

   for i=1:3
       uu             = xx;        % needs to be a matrix for interp2
       uu(:)          = u(i:3:length(u),1);
       data           = interp2(x,y,uu,xxi,yyi);
                        if strcmp(CoordSystem,'Disloc') data=flipud(data); end
       enu(i).data    = data;
       enu(i).x_first = basemap.x_first;
       enu(i).y_first = basemap(1).y_first;
       enu(i).x_step  = basemap(1).x_step ;
       enu(i).y_step  = basemap(1).y_step ;
       enu(i).x_unit  = basemap(1).x_unit ;
       enu(i).Unit    = modelopt.Unit   ;
   end
   %plotdataopt.basemap=basemap;
   %plotdataopt.colorbaropt.Location='OutsideLowerRight'
   %PlotData(enu(3),plotdataopt);
%%
   

   enu(1).DataSet = 'East' ;
   enu(2).DataSet = 'North';
   enu(3).DataSet = 'Up'   ;

   if exist('x_unit') && strcmp(x_unit,'degrees')
       coord=lola2xy(coord,basemap,-1);
   end
