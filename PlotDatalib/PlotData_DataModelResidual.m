function []=PlotData_DataModelResidual(data,opt,rows)
% PlotData_DataModelResidual  -  calls PlotData such that Data, Model and Residual is plotted in 3 rows             
%
%  usage:  PlotData_DataModelResidual(data,opt)
%  TODO: Needs 4th row to plot profile
if ~exist('opt','var')  opt =[]; end
if ~exist('rows','var') rows=4 ; end

   N_SAR     = sum([data(:).SAR]);
   N_GPSonly = sum([data(:).GPSonly]);

   Npx=round(N_SAR/rows)+N_GPSonly;       %subplot windows horizontal
   Npy=rows;                              %subplot windows vertical
   opt.Quakes  = 'off' ;
   opt.gpsfile = 'off' ;
   opt.Coord   = 'off';

   if N_GPSonly && ~N_SAR
      Npy = 1;
   end

   for i=1:length(data)
           tmpopt=opt;
           if i< length(data) tmpopt.colorbaropt.Location='off'; end
           if i==length(data) tmpopt.colorbaropt.Location='OutsideLowerRight'; end
           if strfind(data(i).DataSet,'Topography')  tmpopt.CLim=0; end
           %remove vectors unless we plot Data or Model
           if isempty(strfind(data(i).DataSet,'Data')) && isempty(strfind(data(i).DataSet,'Model')) && isempty(strfind(data(i).DataSet,'Topography'))
                       tmpopt.VectorsBlack=0; tmpopt.VectorsRed=0; 
           end 
           %if isempty(strfind(data(i).DataSet,'Data')) && isempty(strfind(data(i).DataSet,'Model'))  tmpopt.VectorsBlack=0; tmpopt.VectorsRed=0; end 

       eval(sprintf( 'subplot(Npy,Npx,%d)' ,i))

       PlotData(data(i),tmpopt)
   end
