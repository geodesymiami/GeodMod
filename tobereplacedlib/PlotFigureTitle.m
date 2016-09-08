function []=PlotFigureTitle(str)
% PlotFigureTitle      -  plots title on figure             
%     function to plot title on figure   
%     
%     FA, Apr 2004
%
      axes
      th=text(0,0.5,str,'Units','normalized','FontSize',14,'FontName','Courier','HorizontalAlignment','center');
      set(gca,'NextPlot','add','position',[0.5 0.9 0.1 0.10],'vis','off');
      %get(th);


