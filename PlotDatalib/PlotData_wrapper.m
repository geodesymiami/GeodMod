function []=PlotData_wrapper(data,opt)
% PlotData  -  calls PlotData multiple times depending on size of data             
%
%  usage:  PlotData_wrapper(data,opt)
if ~exist('opt','var') opt=[]; end

if length(data)==1
   PlotData(data,opt)
elseif length(data)==2
   subplot(1,2,1)
   PlotData(data(1),opt)
   subplot(1,2,2)
   PlotData(data(2),opt)
else
   N_igrams=length(data);
   Npx=ceil(sqrt(N_igrams));       %subplot windows horizontal
   Npy=ceil((N_igrams/Npx));       %subplot windows vertical
   opt.Coord='off';

   for i=1:N_igrams
       eval(sprintf( 'subplot(Npy,Npx,%d)' ,i))
       PlotData(data(i),opt)
   end
end
