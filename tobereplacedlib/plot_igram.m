function []=plot_igram(igram,varargin)
% plot_igram  -  plots one or a set of interferograms             
%
%  usage:  plot_igram(igram)
%          plot_igram(igram(2:3))
%
%          Input:   igram:   1xN structure array containing N interferograms with:
%
%          optional arguments (given as 'CLim',[-1,1],'LocInd',[2,2],...)
%
%          'CLim'         Color Limits   [ -1,1 ]            [default min and max of data]
%          'CSpaceLim'    window of the interferograms where to compute Color Limits  [default all the image]
%          'Faults'       plot linefile  [ Lllh   (:,2)  array with x,y coordinates ]
%          'Coord'        plot coordinates ['on','off']  [default 'on']
%          'FlipScale'    Flip ColorScale  ['on','off']  [default 'off']
%          'FigName'      Figure Names                   [default 'Interferograms']
%          'ULabel'       Upper Label  'date' or 't'     [default 'date']  
%          'Triangles'    plots triangles                [default none] 
%          'Vectors'      plots vectors                  [default none] 
%          'VertVectors'  plots vertical vectors         [default none] 
%          'VertVectors2' plots 2nd vertical vectors         [default none] 
%          'HorzVectors'  plots horizontal vectors       [default none] 
%          'HorzVectors2' plots 2nd horizontal vectors   [default none] 
%          'Sites'        add GPS site name              [default none] 
%          'Contour'      draw contours                  [default none] 
%          'LocInd'       [50,80]   plots point          [default none] (not tested)
%          'Modulo'     (not tested)  
%
%  Part of the TimeSeries suite
%  FA, March 2005,   


%
% process the arguments
%
if length(varargin) >= 1
if (floor(length(varargin)/2)~=length(varargin)/2) , error ('argument missing') , end
for i=1:2:length(varargin)
    if     strmatch('CLim',varargin{i},'exact')         CLim=varargin{i+1} ;
    elseif strmatch('CSpaceLim',varargin{i},'exact')    CSpaceLim=varargin{i+1}  ;
    elseif strmatch('Faults',varargin{i},'exact')       Faults=varargin{i+1}  ;    Coord='on' ;
    elseif strmatch('Coord',varargin{i},'exact')        Coord=varargin{i+1}  ;
    elseif strmatch('FlipScale',varargin{i},'exact')    FlipScale=varargin{i+1}  ;
    elseif strmatch('Fac',varargin{i},'exact')          Fac=varargin{i+1}  ;
    elseif strmatch('ULabel',varargin{i},'exact')       ULabel=varargin{i+1} ;
    elseif strmatch('FigName',varargin{i},'exact')      FigName=varargin{i+1};
    elseif strmatch('LocInd',varargin{i},'exact')       LocInd=varargin{i+1} ;
    elseif strmatch('Modulo',varargin{i},'exact')       Modulo=varargin{i+1} ;
    elseif strmatch('displot',varargin{i},'exact')      dislocation=varargin{i+1} ;
    elseif strmatch('Triangles',varargin{i},'exact')    Triangles=varargin{i+1} ;
    elseif strmatch('Vectors',varargin{i},'exact')      Vectors=varargin{i+1} ;
    elseif strmatch('VertVectors',varargin{i},'exact')  VertVectors=varargin{i+1} ;
    elseif strmatch('VertVectors2',varargin{i},'exact') VertVectors2=varargin{i+1} ;
    elseif strmatch('HorzVectors',varargin{i},'exact')  HorzVectors=varargin{i+1} ;
    elseif strmatch('HorzVectors2',varargin{i},'exact') HorzVectors2=varargin{i+1} ;
    elseif strmatch('Sites',varargin{i},'exact')        Sites=varargin{i+1} ;
    elseif strmatch('Contour',varargin{i},'exact')      Contour=varargin{i+1} ;
    else
        error('unknown argument: %s',varargin{i});
    end 
end
end
%
% default arguments
%
sc=100000; %scale for quivers
if exist('FigName')==0  FigName='Interferograms'; end
if exist('Coord')==0    Coord='off'; end
if exist('ULabel' )==0  ULabel='date'; end
                % get minimun and maximum displacement rate for plotting if CLim not given
if exist('CLim','var')==0
  if exist('CSpaceLim')
    Xinterv=[CSpaceLim(1):CSpaceLim(2)];Yinterv=[CSpaceLim(3):CSpaceLim(4)]; 
    for i=1:length(igram)
      maxlist(i)=mmax(igram(i).data(Yinterv,Xinterv)) ;
      minlist(i)=mmin(igram(i).data(Yinterv,Xinterv)) ;
    end
  else
    for i=1:length(igram)
      maxlist(i)=mmax(igram(i).data) ;
      minlist(i)=mmin(igram(i).data) ;
    end
  end
   CLim=[min(minlist),max(maxlist)] 
end
if exist('Fac')  
   for i=1:length(igram) igram(i).data=igram(i).data*Fac;  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% start plotting program %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_igrams=length(igram);
Npx=ceil(sqrt(N_igrams));       %subplot windows horizontal
Npy=ceil((N_igrams/Npx));       %subplot windows vertical 

% label vectors for 'Coord'
xl=linspace(igram(1).x_first,igram(1).x_first+igram(1).x_step*size(igram(1).data,2),size(igram(1).data,2));
yl=linspace(igram(1).y_first,igram(1).y_first+igram(1).y_step*size(igram(1).data,1),size(igram(1).data,2));
if exist('dislocation')
   dislocation(6)=igram(1).x_first + dislocation(6)*1000 ;
   dislocation(7)=igram(1).y_first + igram(1).y_step*size(igram(1).data,1) + dislocation(7)*1000 ;
dislocation(1)=dislocation(1)*1000 ;
dislocation(2)=dislocation(2)*1000 ;
dislocation(3)=dislocation(3)*1000 ;
end

%figure('inverthardcopy','off','Color',[1 1 1]);
%figure
%figure('Name',FigName,'Numbertitle','on')
colormap('default');q=colormap(jet);cm=[1 1 1; q(33:61,:);q(6:32,:) ] ;
colormap('default');q=colormap(jet);cm=[q(33:61,:);q(6:32,:) ] ;
q=colormap(jet);jjet=[1 1 1; q(6:31,:);q(32:61,:) ] ;
colormap(cm);colormap(jet)  % use cm for Modulo
for i=1:N_igrams
  titlestr=[];
  eval(sprintf( 'subplot(Npy,Npx,%d)' ,i))   
           if ~isempty(strmatch('date',ULabel)) && isfield(igram,'date1')  titlestr=sprintf('#%d %s-%s',i,igram(i).date1(3:end),igram(i).date2(3:end)) ; end
           if ~isempty(strmatch('t',ULabel))    && isfield(igram,'t1')     titlestr=sprintf('#%d dates:%d-%d',i,igram(i).t1,igram(i).t2) ; end 
           if isfield(igram,'rating') && isfield(igram,'rating') titlestr=strcat(titlestr,sprintf('   r%d',igram(i).rating)); end
  if exist('FlipScale')  igram(i).data=igram(i).data*-1; end
  disp(titlestr)
  if exist('Modulo')
     qqa=mod(igram(i).data,Modulo);
     CLim=[mmin(qqa) mmax(qqa)];
        if exist('Coord')
           imagesc(xl,yl,qqa,CLim)  ;  axis equal ; axis tight ; axis xy;  title(titlestr);
        else
           imagesc(igram(i).data,CLim)  ;  axis equal ; axis tight ; axis off ;  title(titlestr);
        end
     colormap(cm);
  else
        if strmatch(Coord,'on')
           imagesc(xl,yl,igram(i).data,CLim)  ;  axis equal ; axis tight ; axis xy;  title(titlestr);
        else
           imagesc(igram(i).data,CLim)  ;  axis equal ; axis tight ; axis off ;  title(titlestr);
        end
  end
  if exist('Contour','var')  
              qq=medfilt2(igram(i).data,[7 7]);
             [C,h]=contourf(xl,yl,qq,Contour);             ;
             set(h,'ShowText','on','TextStep',get(h,'LevelStep'));
   end
  if exist('dislocation')  hold on ; displot(dislocation) ;  hold off ; end
  if exist('Faults')       hold on ; plot(Faults(:,1),Faults(:,2),'w') ; hold off ; end
  if exist('Triangles')    hold on ; plot(Triangles(:,1),Triangles(:,2),'k^') ; hold off ; end
  if exist('Vectors')      hold on ; quiver(Vectors(:,1),Vectors(:,2),Vectors(:,3),Vectors(:,4),'k') ; hold off ; end
  if exist('HorzVectors')  hold on ; quiver(HorzVectors (:,1),HorzVectors (:,2),HorzVectors (:,3)*sc,HorzVectors (:,5)*0 ,0,'k'); hold off ; end
  if exist('VertVectors')  hold on ; quiver(VertVectors (:,1),VertVectors (:,2),VertVectors (:,3)*0 ,VertVectors (:,5)*sc,0,'r')  ; hold off ; end
  if exist('HorzVectors2') hold on ; quiver(HorzVectors2(:,1),HorzVectors2(:,2)-igram(1).y_step*5 ,HorzVectors2(:,3)*sc,HorzVectors2(:,5)*0,0,'b'); hold off ; end
  if exist('VertVectors2') hold on ; quiver(VertVectors2(:,1)+igram(1).x_step*5,VertVectors2(:,2), VertVectors2(:,3)*0,VertVectors2(:,5)*sc,0,'b')  ; hold off ; end
  if exist('Sites')        hold on ; text  (Sites.xy(1,:),Sites.xy(2,:),Sites.sites(:,:),'FontSize',6,'VerticalAlignment','bottom','Color','k'); hold off ; end


  if exist('LocInd')      hold on ; plot(locind(1),locind(2),'kx') ; hold off ; end
  CBar='on';
  if exist('CBar')  colorbar('horiz') ; end
end

if N_igrams==1 impixelinfo ; end
PlotFigureTitle(FigName)
