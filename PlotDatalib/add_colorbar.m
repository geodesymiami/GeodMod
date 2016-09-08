function []=add_colorbar(opt)
% add_colorbar           -   adds a colorbar to plot  
%
%  TODO; cleaner programming would be:
%        1) call the function as add_colorbar(h,colorbaropt)  where h is the axes handle to which the colorbar is added.
%        2) have colorbaropt default options defined here (in typical default_options style)   
%
%  FA NOV 2006

f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end

   if colsat && length(colsat)>1
      minc = colsat(1);
      maxc = colsat(2);
    else
      minc = min(data(:));
      maxc = max(data(:));
    end

    maxv=max(abs(minc),abs(maxc));
    if  floor(maxv*1000 ) > 0  fac=1000  ; end
    if  floor(maxv*100  ) > 0  fac=100   ; end
    if  floor(maxv*10   ) > 0  fac=10    ; end
    if  floor(maxv*1    ) > 0  fac=1     ; end
    if  floor(maxv*0.1  ) > 0  fac=0.1   ; end
    if  floor(maxv*0.01 ) > 0  fac=0.01  ; end
    if  floor(maxv*0.001) > 0  fac=0.001 ; end
    if  maxc             == 0  fac=1     ; end
    
    xtick_first=floor(abs(minc)*fac)/fac*sign(minc) ;
    xtick_last =floor(abs(maxc)*fac)/fac*sign(maxc) ;
    span_data  = abs(minc)+abs(maxc);
    if span_data==0 span_data=1; end                     % avoids error messages

    image_position=get(gca,'Position') ;

    if strcmp(Location,'InsideLowerLeft')
       xll = image_position(1);
       yll = image_position(2);
    elseif strcmp(Location,'InsideLowerRight')
       xll = image_position(1)+image_position(3)-width-2*x_delta;
       yll = image_position(2);
    elseif strcmp(Location,'OutsideLowerLeft')
       xll = image_position(1);
       yll = 0+0*y_delta;
    elseif strcmp(Location,'OutsideLowerRight')
       xll = image_position(1)+image_position(3)-width-2*x_delta;
       yll = 0+0*y_delta;
    elseif ~ischar(Location)
       xll = Location(1);
       yll = Location(2);
    end

    position_background      = [xll                  yll                         width+2*x_delta height+2*y_delta];
    position_colorbar        = [xll+x_delta          yll+y_delta                 width           height          ];
    position_colorbartitle   = [position_colorbar(1) position_colorbar(2)+height position_colorbar(3:4)          ];

    axes('Position',position_background,'XTick',[],'YTick',[],'XColor',[0 0 0],'YColor',[0 0 0],'Color',[1 1 1],'Box','on')
    %htmp=colorbar;xticks=get(htmp,'XTick');xticklabels=get(htmp,'XTickLabel'); colorbar('off');

    colormap(Cmap);
    cb=colorbar('Location','SouthOutside','Position',position_colorbar) ;
    
             XLim       = get(cb,'XLim') ;
             xticklabels= linspace (xtick_first,xtick_last,label_number) ;
             xticks     = (xticklabels-minc)/span_data * (XLim(2) - XLim(1))+XLim(1) ;
             %set(cb,'XLim',[minc,maxc])
             set(cb,'XTick',xticks,'XTickLabel',xticklabels,'FontSize',8)
    annotation('textbox',position_colorbartitle,'EdgeColor','none','HorizontalAlignment','Center','String',Title,'FontSize',8);
