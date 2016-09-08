function []=PlotData(data,opt)
% PlotData  -  plots one or a set of interferograms             
%
%  usage:  PlotData(data,opt)
%
%          Input:   data:   1x1 structure array containing N interferograms with:
%
%          optional arguments (given as 'CLim',[-1,1],'LocInd',[2,2],...)
%
%          'Shade'          Input Shaded relief (.jpg only)                             [default 'off']
%          'ShadeOnly'      Plot shaded relief without DEM                              [default 'off']
%          'CLim'           Color Limits   ([-1 1],'Centered')                          [default min and max of data]
%          'Cmap'           'dismph, 'blueredN', 'jet'                                  [default  modified jet]
%          'CSpaceLim'      window of the interferograms where to compute Color Limits  [default all the image]
%          'Faults'         plot linefile  [ Lllh   (:,2)  array with x,y coordinates ]
%          'Coord'          plot coordinates ['on','off']                               [default 'on']
%          'dopixval'       Get value using pixval
%          'FlipColorScale' Flip ColorScale  ['on','off']  (note that the default       [default 'off']
%                           is flipud(jet) so that uplift (negative) has warm colors)    
%          'ULabel'         Title of Plot                                               [default 'Date']  
%          'Triangles'      plots triangles                                             [default 'off'] 
%          'VectorsBlack'   plots vectors                                               [default 'off'] 
%          'VectorsRed'    plots vectors                                               [default 'off'] 
%          'VertVectors'    plots vertical vectors                                      [default 'off'] 
%          'VertVectors2'   plots 2nd vertical vectors                                  [default 'off'] 
%          'HorzVectors'    plots horizontal vectors                                    [default 'off'] 
%          'HorzVectors2'   plots 2nd horizontal vectors                                [default 'off'] 
%          'Sites'          add GPS site name                                           [default 'off'] 
%          'Quakes'         add Quakes                                                  [default 'off'] 
%          'Focals'         add Focals                                                  [default 'off']
%          'PlotContour'    struct('values',[],'color','k.','label','on' draw contours  [default 'off'] 
%          'LocInd'         [50,80]   plots point                                       [default 'off'] (not tested)
%          'Symbols'        user defined plot symbols [e.g. symbols.Faults={'k'}]       [default 'off'] 
%          'Fringe'         color cycle spacing (1: calculate from wavelength,time)      [default 'off'] 
%          'ShadeFac'       factor for ration Data/Shade                                [default 0]
%          'PlotUnit'       'Radian','m/yr','cm/yr','mm/yr'                             [default data.Unit]
%                           data are scaled if necessary   
%          'colorbaropt'    structure with colorbaroptions 
%                           'Location' [xll yll],'InsideLowerLeft' ,'InsideLowerRight', [default 'OutsideLowerRight']
%                                                'OutsideLowerLeft','OutsideLowerRight'
%                           'Title'                                                     [default 'LOS velocity']
%                           'label_number' number of xlabels                            [default  5]
%                           'width'   in figure coordinates                             [default  0.16]
%                           'height'  in figure coordinates                             [default  0.03]
%                           'x_delta' extension of axes for white background            [default  0.03]
%                           'y_delta' extension of axes for white background            [default  0.04]
%
%  Part of the TimeSeries suite
%  V1.0  October  2005:  Falk Amelung and Noel Gourmelen, October 2005   
%  V1.01 November 2005:  added Quakes and Symbol structure functionality   (FA)
%  V1.02 November 2005:  added ShadeFac functionality   (NG)'/RAID6/insar_lab/testdata_geodmod/Wells/Data/EnvD2'
%  V1.03 September2006:  Addjusted for MakeDataModel suite  (FA)
%  V1.04 November 2006:  added colorbar option  (FA)

%  TODO: Assign Quakes.lola or Quakes.xy to Quakes.coord depending of x_unit=degree or km
%  Write a function subset_vectordata to extract data for area and depth of
%  interest, e.g.:
%  Quakes=extract_subset(Quakes,basemap,x_unit,'Quakes',dep_ranges)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultopt=struct(                                        ... 
        'PlotType'       ,        '2D'       ,            ...
        'CLim'           ,        'off'      ,            ...
        'CSpaceLim'      ,        'off'      ,            ...
        'ColorTable'     ,        'off'      ,            ...
        'Cmap'           ,        'off'      ,            ...
        'CBar'           ,        'on'       ,            ...
        'dopixval'       ,        'off'      ,            ...
        'PlotContour'    ,        'off'      ,            ...
        'FlipColorScale' ,        'off'      ,            ...
        'Coord'          ,        'on'       ,            ...
        'ShadeOnly'      ,        'off'      ,            ...
        'Dem'            ,        'off'      ,            ...
        'Dislocation'    ,        'off'      ,            ...               
        'modelopt'       ,        'off'      ,            ...               
        'Fac'            ,        'off'      ,            ...
        'ULabel'         ,        'off'      ,            ...
        'DataName'       ,        'off'      ,            ...
        'LocInd'         ,        'off'      ,            ...
        'Fringe'         ,        'off'      ,            ...
        'Displot'        ,        'off'      ,            ...
        'Sites'          ,        'off'      ,            ...
        'Faults'         ,        'off'      ,            ...
        'Roads'          ,        'off'      ,            ...
        'Quakes'         ,        'off'      ,            ...
        'Profile'        ,        'off'      ,            ...
		'Focals'         ,        'off'      ,            ...
        'GPSdata'        ,        'off'      ,            ...
        'GPSpred'        ,        'off'      ,            ...
        'Triangles'      ,        'off'      ,            ...
        'HorzVectors'    ,        'off'      ,            ...
        'HorzVectors2'   ,        'off'      ,            ...
        'Symbols'        ,        'off'      ,            ...
        'VectorsBlack'   ,        'off'      ,            ...
        'VectorsRed'     ,        'off'      ,            ...
        'VertVectors'    ,        'off'      ,            ...
        'ShadeFac'       ,         0.2       ,            ...
        'PlotUnit'       ,        'off'      ,            ...
        'x_unit'         ,        'degrees'  ,            ...
        'Surface3D'      ,        'off'      ,            ...
        'z_offset'       ,        'off'      ,            ...
        'VertExaggeration',          1       ,            ...
        'DownSampleFac'  ,        'off'      ,            ...
        'viewdir'        ,        'off'      ,            ...
        'cam_position'   ,        'off'      ,            ...
        'alpha_val'      ,           1       ,            ...
        'TopoMedfilt'    ,          11       ,            ...
        'vec2bitmap'     ,           0       ,            ...
        'marker1D'       ,        'off'      ,            ...
        'VertVectors2'   ,        'off'      )            ;
defaultopt.colorbaropt=struct(                            ... 
        'Location'       ,        'InsideLowerRight' ,     ...
        'Title'          ,        'off'             ,     ...
        'width'          ,        0.16              ,     ...
        'height'         ,        0.03              ,     ...
        'x_delta'        ,        0.03              ,     ...
        'y_delta'        ,        0.04              ,     ...
        'label_number'   ,         5         )            ;
defaultopt.googleearthopt=struct(                         ... 
        'DoIt'           ,        'off' ,                  ...
        'option1'        ,        'off'    ,     ...
        'option2'        ,         0         )            ;
defaultopt.distribopt    =struct(                         ...
        'PlotThresh'     ,        5          )            ;
if ~exist('opt','var')  opt=[]; end
[opt]=process_defaultoptions(opt,defaultopt);  %display(opt)
%[opt.googleearthopt]=process_defaultoptions(opt.googleearthopt,defaultoptr.googleearthopt); display(googleearthopt)
if ~isfield(data,'data') data.data=opt.basemap.data; data.amp=ones(size(data.data)) ; end   % this is useful if GPS data are given
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Some default values  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sc=100000; %scale for quivers
if ~ ShadeFac        ShadeFac=0.0 ;               end
if   ShadeOnly       colorbaropt.Location=false ; end   
[symbols]=process_symboloptions(opt,Symbols) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% label vectors for 'Coord' Quakes.xy(:,1)
%
 x_l=linspace(data(1).x_first,data(1).x_first+data(1).x_step*size(data(1).data,2),size(data(1).data,2));
 y_l=linspace(data(1).y_first,data(1).y_first+data(1).y_step*size(data(1).data,1),size(data(1).data,1));

 if strcmp(x_unit,'km')
     x_l = lola2xy([x_l'   x_l'*0]',basemap,1); x_l = x_l(1,:);
     y_l = lola2xy([y_l'*0 y_l'  ]',basemap,1); y_l = y_l(2,:);
 end

 if strcmp(x_unit,'degrees') || strcmp(x_unit,'degres')
     what = 'lola';
 elseif strcmp(x_unit,'km')
     what='xy';
 end
 
 if isstruct(Quakes)            Quakes.X      = Quakes.(what)(:,1);      Quakes.Y      = Quakes.(what)(:,2);     end
 if isstruct(Focals)            Focals.X      = Focals.(what)(:,1);      Focals.Y      = Focals.(what)(:,2);     end
 if isstruct(Faults)            Faults.X      = Faults.(what)(:,1);      Faults.Y      = Faults.(what)(:,2);     end
 if isstruct(Roads)             Roads.X       = Roads.(what)(:,1);       Roads.Y       = Roads.(what)(:,2);      end
 if isstruct(Profile)           Profile.X     = Profile.(what)(:,1);     Profile.Y     = Profile.(what)(:,2);    end
 if isstruct(GPSdata)           GPSdata.X     = GPSdata.(what)(1,:);     GPSdata.Y     = GPSdata.(what)(2,:);    end
 if isstruct(GPSpred)           GPSpred.X     = GPSpred.(what)(1,:);     GPSpred.Y     = GPSpred.(what)(2,:);    end 
 if iscell(marker1D)
    for i=1:length(marker1D)    marker1D{i}.X = marker1D{i}.(what)(1,:); marker1D{i}.Y = marker1D{i}.(what)(2,:);end
 end
 %end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Prepare for plotting %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~PlotUnit
   PlotUnit=data.Unit;
end
if  strcmp(PlotUnit,'radian')  ||  strcmp(PlotUnit,'m/yr') || strcmp(PlotUnit,'Meters')  Fac = 1    ;
    elseif strcmp(PlotUnit,'unity')                                                      Fac = 1    ;
    elseif strcmp(PlotUnit,'m')                                                          Fac = 1    ;
    elseif strcmp(PlotUnit,'cm/yr')                                                      Fac = 100  ;
    elseif strcmp(PlotUnit,'mm/yr')                                                      Fac = 1000 ;
    else error('user error: PlotUnit %s not defined -- exiting',PlotUnit)
end
data.data=data(1).data*Fac;    

[titlestr,colorbar_title]=GenerateNamesFromDataStructure(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Plot Data %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Fringe
    load dismph
    colormap(dismph);
    if Fringe==1                                     % calulate default fringe cycle from wavelength 
        Fringe = CalcDefaultFringe(data,PlotUnit);   % (and from TotalTime if PlotUnit m/yr)
    end
    ph     = mod(data.data,Fringe);
    CLim   = false ;
else
    ph=data.data;
end

if CLim
  if strcmp(CLim,'Centered')
     maxabs=max(abs(mmax(data(1).data)),abs(mmin(data(1).data))) ;
     CLim=[-maxabs maxabs];
  end
     
  ph(find(ph<CLim(1)))=CLim(1);
  ph(find(ph>CLim(2)))=CLim(2);
  if isempty(find(ph<CLim(1)))
     ph(1,2)=CLim(1);
  end
  if isempty(find(ph>CLim(2)))
     ph(1,1)=CLim(2);
  end
  %whos ph
end 

if exist('sun_direction')==0
    sun_direction=[60 45];
end

if length(Cmap)-1
  load(Cmap);
  Cmap=eval(Cmap);
else
  Cmap=jet;
end
if ~FlipColorScale
  %Cmap=flipud(Cmap);
end

if strcmp(data.DataSet,'Topography')   % sbaker 10/31/2008, then again on 12/07/2008 because it was lost
    Cmap = demcmap(Dem);               % Use matlabs built-in colormap for topgraphy/bathymetry
end

  hold all;
  if isstruct(PlotContour)
      figure;
      set(gcf, 'BackingStore', 'off');
      set(gcf, 'renderer', 'OpenGL');
      axes('position', [0,0,1,1]);
      
      %%% CHANGES THESE TO GET GOOD CONTOURS FOR POS AND NEG VALUES
%       contourf(x_l,y_l,ph,'LevelList',opt.PlotContour.values,'LineStyle','none')
%       load '/RAID6/sbaker/models/Kilauea2005-2007/Contours/contourPos.mat';
%       load '/RAID6/sbaker/models/Kilauea2005-2007/Contours/jetPos.mat';
%       Cmap=contourPos; Cmap=jetPos;
      contourf(x_l,y_l,ph*-1,'LevelList',opt.PlotContour.values,'LineStyle','none')
      load '/RAID6/sbaker/models/Kilauea2005-2007/Contours/contourNeg.mat';
      load '/RAID6/sbaker/models/Kilauea2005-2007/Contours/jetNeg.mat';
      Cmap=contourNeg; Cmap=jetNeg;
      
      caxis(CLim);
      %Cmap=colormap(blueredN);
      colormap(Cmap);
      axis tight;
      set(gca,'LooseInset',get(gca,'TightInset'))
      axis off;
      tmpImage=imcapture(gcf,'img',[size(ph,1) size(ph,2)]);
      close gcf;

      ph2=rgb2ind(tmpImage,Cmap);
      ph2=double(ph2);
      for ni=0:length(Cmap)
      ph2(find(ph2==ni)) = CLim(1)+((CLim(2)-CLim(1))/length(Cmap))*ni;
      end
      ph2(find(isnan(ph)==1))=NaN;
      ph=ph2;clear ph2;
% =======
%      qq=medfilt2(ph,[7 7]);
%      wwidth=size(PlotContour.values,2);
%      valuesC=PlotContour.values;
%      Contours(1)=valuesC(1)-(valuesC(2)-valuesC(1));Contours(2:wwidth+1)=valuesC;Contours(wwidth+2)=valuesC(wwidth)+(valuesC(wwidth)-valuesC(wwidth-1));
%      ph2=zeros(size(ph));
%      for i=1:size(Contours,2)-1
%          ph2(find(ph<Contours(i+1) & ph>=Contours(i)))=Contours(i);
%      end
%      ph2(find(isnan(ph)==1))=NaN;
%      ph=ph2;clear ph2;
% >>>>>>> 1.16
  end

  if isstruct(PlotContour)
    if isfield(PlotContour,'color')
       [C,h]=contour(xl,yl,ph,Contours); 
       plot(C(1,:),C(2,:),PlotContour.color,'MarkerSize',1); 
       if isfield(PlotContour,'label')
           clabel(C,h);
       end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%% Add symbols, faults, etc. %%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isstruct(Quakes)
       for j=length(symbols.Quakes.mag_ranges)-1:-1:1
	       ind=find(Quakes.mag >= symbols.Quakes.mag_ranges(j) & Quakes.mag < symbols.Quakes.mag_ranges(j+1)) ;
           tmpQuakes.X     = Quakes.X(ind);
           tmpQuakes.Y     = Quakes.Y(ind);
           tmpQuakes.depth = Quakes.depth(ind);
           tmpQuakes.mag   = Quakes.mag(ind);
           for i=length(symbols.Quakes.dep_ranges)-1:-1:1
	           ind=find(tmpQuakes.depth >= symbols.Quakes.dep_ranges(i) & tmpQuakes.depth < symbols.Quakes.dep_ranges(i+1)) ;
                       symb       = {symbols.Quakes.dep_symbols{i},symbols.Quakes.mag_symbols{j}{:} } ;
                       symb       = {symb{:} 'MarkerFaceColor',sprintf('%.1s',symb{1})};
                       plot(tmpQuakes.X,tmpQuakes.Y,symb{:});  
           end
       end
    end

    if isstruct(Focals)
           tmpFocals.X     = Focals.X;
           tmpFocals.Y     = Focals.Y;
           tmpFocals.depth = Focals.depth;
           tmpFocals.mag   = Focals.mag;
		   tmpFocals.strike= Focals.strike;
	       tmpFocals.dip   = Focals.dip;
		   tmpFocals.rake  = Focals.rake;
		   tmpFocals.size  = Focals.size;
           hold on
           for h=1:size(tmpFocals.X,1);
		   beachball(tmpFocals.strike(h),tmpFocals.dip(h),tmpFocals.rake(h),tmpFocals.X(h),tmpFocals.Y(h),tmpFocals.size(h));
           end
    end

    % note: use Faults(1) etc because some values are NaN. No NaN's should be fed into this.
    if isstruct(modelopt) plot_model_parameters(modelopt,x_unit,PlotType)                                                                    ;end
    if isstruct(Profile)  plot(Profile.X,Profile.Y,'--k')                                                                                     ;end
    if Triangles(1)    plot(Triangles(:,1),Triangles(:,2),symbols.Triangles{:})                                                              ;end
    if isstruct(GPSdata)      quiver(GPSdata.X,GPSdata.Y,GPSdata.enu(1,:),GPSdata.enu(2,:),symbols.VectorsBlack{:})                          ;end   
    if isstruct(GPSpred)      quiver(GPSpred.X,GPSpred.Y,GPSpred.enu(1,:),GPSpred.enu(2,:),symbols.VectorsRed{:})                            ;end    
   
    GPSvert=1; % FA 12/08  needs to be carried in through options. For now modify manually
    if GPSvert
       if isstruct(GPSdata)      quiver(GPSdata.X,GPSdata.Y,zeros(1,length(GPSdata.X)),GPSdata.enu(3,:),symbols.VectorsBlack{:})              ;end 
       if isstruct(GPSdata)      quiver(GPSpred.X,GPSpred.Y,zeros(1,length(GPSpred.X)),GPSpred.enu(3,:),symbols.VectorsRed{:})                ;end
    end
    % FA 12/08 the following lines for Vectors are probably obsolete (later
   % FA 12/08  above comment not true. VectorsBlack is used in
   % PlotTheModel(vertical field plus horizontal vectors). This should be
   % converted to structur similar to GPS so that it works for xy also
    if sum(VectorsBlack(:,1)) quiver(VectorsBlack(:,1),VectorsBlack(:,2),VectorsBlack(:,3),VectorsBlack(:,4),symbols.VectorsBlack{:})        ;end
    %if sum(VectorsRed(:,1))   quiver(VectorsRed  (:,1),VectorsRed  (:,2),VectorsRed  (:,3),VectorsRed  (:,4),symbols.VectorsRed  {:})        ;end
    if sum(HorzVectors(:,1))  quiver(HorzVectors (:,1),HorzVectors (:,2),HorzVectors (:,3)*sc,HorzVectors (:,5)*0 ,0,symbols.HorzVectors{:}) ;end
    if sum(VertVectors(:,1))  quiver(VertVectors (:,1),VertVectors (:,2),VertVectors (:,3)*0 ,VertVectors (:,5)*sc,0,symbols.VertVectors{:}) ;end
    if sum(HorzVectors2(:,1)) quiver(HorzVectors2(:,1),HorzVectors2(:,2)-data(1).y_step*5 ,HorzVectors2(:,3)*sc,HorzVectors2(:,5)*0,0,symbols.HorzVectors2{:}) ;end
    if sum(VertVectors2(:,1)) quiver(VertVectors2(:,1)+data(1).x_step*5,VertVectors2(:,2), VertVectors2(:,3)*0,VertVectors2(:,5)*sc,0,symbols.VertVectors2{:}) ;end
    if isstruct(Sites) text  (Sites.xy(1,:),Sites.xy(2,:),Sites.sites(:,:),'FontSize',6,'VerticalAlignment','bottom','Color',symbols.Sites{:})           ;end
    if LocInd          plot(locind(1),locind(2),symbols.LocInd{:});                                                                end

    % TODO: shadefile should be read in in Prepare. Did not work because preocess_defaultoptions gave an error related to nested structures (checking for 'off' and

    %if shapefile        
     %  faults=shaperead(shapefile);
      % plot([faults.X],[faults.Y],symbols.Faults{:}); 
    %end
    if isstruct(Faults) plot(Faults.X,Faults.Y,symbols.Faults{:}); end
    if isstruct(Roads)  plot(Roads.X, Roads.Y, symbols.Roads{:}); end
    
    h=gca;
    set(h,'XLim',[min(x_l),max(x_l)]);
    set(h,'YLim',[min(y_l),max(y_l)]);

    if ~Coord set(h,'xticklabel','','yticklabel',''); end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Generate rgb image 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ShadeOnly
   [rgbim]=basemap.shade;
else
   [rgbim]=PlotRGB_Demjpg_N(ph,basemap.shade,ShadeFac,Cmap,1.5,[-75 75],1,colorbaropt,CLim);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   daspect([basemap.y_posting basemap.x_posting 1]);                  % 9/2008: last number relates to vert exag.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now generate 2D or 3D image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch PlotType
case {'2D', 'map'}
    
   h1=imagesc(x_l,y_l,rgbim);    uistack(h1,'bottom');
   title(titlestr)

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % add colorbar
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if colorbaropt.Location
      colorbaropt.colsat   = false ;                 % Colorscale saturation is currently not build into code but easy to implement
      colorbaropt.data     = ph ;
      colorbaropt.Cmap     = Cmap;
      if ~colorbaropt.Title   colorbaropt.Title = colorbar_title ; end
      add_colorbar(colorbaropt);
   end

case{'3D'}
    %
    %  3D image
    %
    %vec2bitmap=1; 
    if vec2bitmap
        bitmap=imcapture(h,'img',[size(rgbim,1) size(rgbim,2)]);
        ind = find(sum(bitmap,3)~=3*255); len=size(bitmap,1)*size(bitmap,2);
        qind = [ind; len+ind; 2*len+ind];
        rgbim(qind)=bitmap(qind); %imagesc(rgbim);
     end
    
   
    if TopoMedfilt
        tmp=medfilt2(basemap.data,[TopoMedfilt TopoMedfilt]);
        tmp(find(tmp==0)) = basemap.data(find(tmp==0));
        basemap.data      = tmp;  
    end
    if DownsampleFac
        [basemap.data,rgbim]  = resampledata(DownsampleFac,basemap.data,rgbim);
        [x_l,y_l]             = resampledata(DownsampleFac,x_l,y_l); 
    end
    z = (basemap.data-min(min(basemap.data)))/1000*VertExaggeration;
    z = z - min(min(z));                                           % using min instead of mean so that shallow quakes plot under the ground
    
    delete(h); hold all;
    h=surface(x_l,y_l,z,double(rgbim)/255); shading flat; 
    alpha(alpha_val);
    
    if isstruct(modelopt)
        [s] = plot_model_parameters(modelopt,x_unit,PlotType);      %slip s returned for colorbar
    end
    
    grid on;
    %camproj('perspective');
    xlabel('East [km]'); ylabel('North [km]'); zlabel('Depth [km]');
    if modelopt.N_disloc || modelopt.N_multidisloc  set(gca,'CameraTarget',[modelopt.par.xy(6) modelopt.par.xy(7),-modelopt.par.xy(3)]); end

    if cam_position(1) 
        set(gca,'CameraPosition',cam_position);
    else 
        view(viewdir);
    end

    if isstruct(Quakes)
        symb={'b' 'r'};
        ind=find(Quakes.depth >= symbols.Quakes.dep_ranges(1) & Quakes.depth <= symbols.Quakes.dep_ranges(end)) ; %remove quakes out of depth ram
        Quakes.X     = Quakes.X(ind);
        Quakes.Y     = Quakes.Y(ind);
        Quakes.depth = Quakes.depth(ind);
        Quakes.mag   = Quakes.mag(ind);
        for j=length(symbols.Quakes.mag_ranges)-1:-1:1
	       ind=find(Quakes.mag >= symbols.Quakes.mag_ranges(j) & Quakes.mag < symbols.Quakes.mag_ranges(j+1)) ;
           scatter3(Quakes.X(ind),Quakes.Y(ind),-Quakes.depth(ind),symbols.Quakes.mag_symbols{j}{2}*5,symb{j},'filled'); 
        end
    end
   
    if modelopt.N_disloc
        load slipcol;
        cmap = colormap( [ones(round(64*distribopt.PlotThresh/100),3) ; slipcol] );        %PlotThresh is given in percent
        slipcolorbaropt          = colorbaropt;
        slipcolorbaropt.Location = 'OutsideLowerLeft';
        slipcolorbaropt.colsat   = false ;                 % Colorscale saturation is currently not build into code but easy to implement
        slipcolorbaropt.data     = s;
        slipcolorbaropt.Cmap     = cmap;
        slipcolorbaropt.Title    = ['slip ''[m]'] ;
        add_colorbar(slipcolorbaropt);
    end

case{'SampledData'}
     close;
     figure;imagesc(rgbim); 
     hold on;
     set(0,'DefaultPatchLineStyle','-')
     patch(data.cx,data.cy,data.datavec,'FaceAlpha',0.5);
     caxis(CLim);
     title(['Sampled ' titlestr ])
     hold off;
     set(0,'DefaultPatchLineStyle','none')
     daspect([basemap.y_posting basemap.x_posting 1]); 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % add colorbar
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   colorbaropt.Location ='OutsideLowerRight';
   if colorbaropt.Location
      colorbaropt.colsat   = false ;                 % Colorscale saturation is currently not build into code but easy to implement
      colorbaropt.data     = ph ;
      colorbaropt.Cmap     = Cmap;
      if ~colorbaropt.Title  colorbaropt.Title = colorbar_title ; end
      add_colorbar(colorbaropt);
   end

     %patch(data.cx+x(1),data.cy+Y(1),data.datavec,'FaceAlpha',0.5);   %FA 1/2011 Sjonni had the X(1) ofset in his code. Not sure what is the purpose
case{'1D' 'CrossSection'}
      if ~strcmp(x_unit,'km') errordlg('PlotData 1D works only for x_unit km');error('user error -- exiting'); end % need x_unit = 'km';

      cla reset
      hold all;

      if isstruct(Profile)
           [xx,yy]       = meshgrid(x_l,y_l);
           xp            = linspace(Profile.xy(1,1),Profile.xy(2,1),200) ;
           yp            = linspace(Profile.xy(1,2),Profile.xy(2,2),200) ;
           profile_dist  = sqrt((xp-xp(1)).^2 + (yp-yp(1)).^2);


           d             = interp2(xx,yy,data.data,xp,yp,'nearest');
           hgt           = interp2(xx,yy,(basemap.data-min(min(basemap.data)))/1000*VertExaggeration,xp,yp,'nearest');
           hgt           = hgt - min(hgt);            % using min instead of mean so that shallow quakes plot under the ground

           [ax]=plotyy(profile_dist,hgt,profile_dist,d) ;
           
           if isstruct(Quakes)
           Quakes.profile_dist = Profile.unitvec*[Quakes.X'-Profile.xy(1,1); Quakes.Y'-Profile.xy(1,2)];            
           
           symb={'b' 'r'};
           ind                 = find(Quakes.depth >= symbols.Quakes.dep_ranges(1) & Quakes.depth <= symbols.Quakes.dep_ranges(end)) ; %remove quakes out of depth ram
           Quakes.profile_dist = Quakes.profile_dist(ind);
           Quakes.depth        = Quakes.depth(ind);
           Quakes.mag          = Quakes.mag(ind);
           for j=length(symbols.Quakes.mag_ranges)-1:-1:1
	          ind=find(Quakes.mag >= symbols.Quakes.mag_ranges(j) & Quakes.mag < symbols.Quakes.mag_ranges(j+1)) ;
              plot(Quakes.profile_dist(ind),-Quakes.depth(ind),[symb{j} 'o'],'MarkerSize',symbols.Quakes.mag_symbols{j}{2},'MarkerFaceColor',symb{j});
           end 
           end
           
           plot_model_parameters(modelopt,x_unit,Profile,PlotType);
           
           %axis equal;
           set(ax(1),'YLimMode','auto', 'YTickMode','auto');
           set(ax(1),'XLim',[0 max(profile_dist)]);
           set(ax(2),'XLim',[0 max(profile_dist)]);
           
           axis equal;
           %plot(ax(1),[0 profile_dist(end)],[0 0],'.k')
           
           ylim1     = get(ax(1),'YLim');
           ylim2     = get(ax(2),'YLim');
           ylimratio = ylim1./ylim2;
           
           if  any([ isinf(ylimratio) isnan(ylimratio)])    % FA 12/2008: Got an error here for KilaueaIntrusion.min
               logmessage('###ATTENTION:### unexplained ylimratio=inf,  skipping ylim adjustment, need to investigae');
           else
               set(ax(2),'YLim',ylim1/ylimratio(2));
           end
           
           xlabel(ax(1),['Distance[' x_unit ']']); 
           ylabel(ax(1),['Depth [km]']); 
           ylabel(ax(2),['Vertical motion [' data(1).Unit ']']); 
           
           if iscell(marker1D)
              for i=1:length(marker1D)
                  marker1D{i}.profile_dist = Profile.unitvec*[marker1D{i}.X-Profile.xy(1,1); marker1D{i}.Y-Profile.xy(1,2)];
                  plot(marker1D{i}.profile_dist,0.1,[marker1D{i}.symb]);
              end
           end
          
      end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Google Earth : save kmz file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if googleearthopt.DoIt                   
      %make_kmz(data,rgbim);
end
%  TODO; cleaner programming would be to call PlotData from logplot such
%  that the kml_string and bitmap are returned. The files are then saved in logplot (with the same name as the pdf) and compressed from
%  logplot.
%   
%  [kml_string,kml_bitmap]=PlotData(...)                % in logplot
%  [kml_string]=make_kml_file(googleearthopt)            % in PlotData  
%  [kml_bitmap]=rgbim                                    % in PlotData
%
%   Else, if the contents of the kmz file can be generated directly in matlab (on the fly compression in memory), then 
%  [kmz_data]=PlotData(...)                % in logplot
%  and just save the kmz_data as a binary file
%  TODO: googleearthopt should decide whether vectors/faults are plotted as a bitmap, etc.
return
strcmp(x_unit,'km')

