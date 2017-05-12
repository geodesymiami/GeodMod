function [data_mod,predout,enu,dataset]=PlotModel(par,dataset,inverseopt,plotmodelopt)
%PlotModel   -  Plot interferogram, models and residual
%
%usage: PlotModel(mhat,objfunc,dataset,inverseopt,objfuncopt,plotmodelopt)
%
% Input:   dataset       structure containing input data
%          objfunc       objective function
%          inverseopt    structure with options for inverse modelling,
%                        including objfuncopt
%          plotmodelopt  structure with options for plotting 
%
% fields of plotmodelopt                       
% Plot             'off' for calculation of rms only       [default 'on']
%                  '2D'  2D field for data and model                  
%                  'GPS' GPS data and model                
% RemoveRamps      'on' for removal of phase ramps for plot[default 'off']
% SampleMethod         'Grid' or 'Quadtree'                          [default 'Grid']
% Vectors          'on'                                    [default
%
% TODO: JUNE 2007: NEED TO ELIMINATE readfrom_dataset_structure 
% USE dataset.data_mod, dataset.predvec AND GENERATE dataset.pred
% USING FORWARDMODEL AND INTERP2 
% June 2007. Already have changed GenericObjectiveFunction and ModifyDatasetLin
%            so that dataset.data_mod and dataset.predvec are returned
%
f=fieldnames(plotmodelopt); for i=1:length(f) eval([char(f{i}) '= plotmodelopt.(f{i}) ;' ]) ; end

ProfileFillWholes = false;

objfuncopt = inverseopt.objfuncopt;
objfunc    = inverseopt.objfuncopt.objfunc ;
modelopt   = inverseopt.objfuncopt.modelopt ;
%Profile    = plotmodelopt.plotdataopt.Profile;

if isfield(plotdataopt,'Faults')   Faults.X =plotdataopt.Faults.xy(:,1);  Faults.Y= plotdataopt.Faults.xy(:,2);  else Faults =false; end
if isfield(plotdataopt,'Profile') &&  isfield(plotdataopt.Profile,'xy') clear Profile; Profile.X=plotdataopt.Profile.xy(:,1); Profile.Y= plotdataopt.Profile.xy(:,2); else Profile=false; end
FontSize   = 7;


x_unit     = 'km';           %TODO: PlotModel should accept modelopt as argument and x_unit shoul determine whether to plot ll or km
                             %TODO: instead of displotmulti plot_model_parameter should be used.

plotstr=GenerateSummary(modelopt,dataset,inverseopt);

sqrflag = true;
[resi,pred,u,rms,npar,weights,rms_unitsig,mlin,dataset] = feval(objfunc,par,dataset,objfuncopt,sqrflag);                     % FA Feb 2008. Needed to make it work for ProjectProfile
par     = npar;

if strcmp(dataset(1).CoordSystem,'ProjectProfile')
   data=dataset(:).datavec;  
   disp('NOEL: plot data and pred ! ');
   plotmodelopt.y_unit = dataset.Unit ;
   logplot('PlotProfile',out_name,dataset.coord(1,:),[data;pred'],plotmodelopt);
   whos data pred
   data_mod=[];predout=[];enu=[]; % needed because data_mod is returned
   return
end

[SAR_1,SAR_2,SAR_3,SAR_4,SAR_5,GPShorz,GPSvert]               = deal(false);
[radarlook_1,radarlook_2,radarlook_3,radarlook_4,radarlook_5] = deal([]);  %Sep 9 2005, commented out
[data_1,data_2,data_3,data_4,data_5]                          = deal([]);
[pred_1,pred_2,pred_3,pred_4,pred_5]                          = deal([]);

readfrom_dataset_structure ;
if SAR_1 data_1 = flipud(data_mod_1);  amp_1 = flipud(amp_1); end
if SAR_2 data_2 = flipud(data_mod_2);  amp_2 = flipud(amp_2); end
if SAR_3 data_3 = flipud(data_mod_3);  amp_3 = flipud(amp_3); end
if SAR_4 data_4 = flipud(data_mod_4);  amp_4 = flipud(amp_4); end
if SAR_5 data_5 = flipud(data_mod_5);  amp_5 = flipud(amp_5); end

pred = [dataset(:).predvec];

if GPShorz || GPSvert
   [dGPS,mGPS] = datasetstruct2GPSstruct(dataset);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove FactorNonLin model parameters from par for plotting (note: introduce objfuncopt.partype)
if strcmp('SAR',objfuncopt.FactorNonLin) || strcmp('GPS',objfuncopt.FactorNonLin) par=par(1:end-1) ; end
if strcmp('SARmul',objfuncopt.FactorNonLin) iSAR=SAR_1+SAR_2+SAR_3+SAR_4+SAR_5; ifac= iSAR-1 + (GPShorz || GPSvert); par=par(1:end-ifac); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_l=([1:size(inverseopt.plotdataopt.basemap.data,2)] - 1)*inverseopt.plotdataopt.basemap.x_posting;   %axis label vector
y_l=([1:size(inverseopt.plotdataopt.basemap.data,1)] - 1)*inverseopt.plotdataopt.basemap.y_posting;   %axis label vector
[xx,yy]=meshgrid(x_l,y_l);

% PLOT
ClimQt=[-0.05 0.02]          ;   ClimQt_res=ClimQt*0.5;
ClimQt=[-0.1 0.1]            ;   ClimQt_res=ClimQt*0.5;
%modulo=0.2                   ;   modulo_res=0.2 ;
%modulo=0.03                  ;   modulo_res=0.03 ;
%modulo=0.50                  ;   modulo_res=0.50;
%modulo=0.028                  ;   modulo_res=0.028 ;
%if strcmp(dataset(1).DataSet,'JersD') |   strfind(dataset(1).DataSet,'Alos')  modulo=0.12 ;  modulo_res=0.12 ; end

% FA 6/2010: modulo is now calculated from wavelength unless specified in inverseopt.plotmodelopt.modulo
% Todo: Specify plotmodelopt independeently of inverseopt but this requires changes in geodmod.m
%
PlotUnit='m';         %FA 6/2010 (don't know how to deal with case of rate. I guuess the unit should be carried into PLotModel

if dataset(1).SAR
   Fringe = CalcDefaultFringe(dataset(1),PlotUnit);
else
   Fringe = 0.1;   % plot deformation in 10 cm fringes if only GPS is given
end

if ~modulo      modulo = Fringe;     end
if ~modulo_res  modulo_res = Fringe; end

%modulo = modulo * 4;       % FA 12/2016  test to improve CSK plots

ClimGrid=[-modulo/56 modulo] ;   ClimGrid_res=[-modulo_res/56 modulo_res];
%modulo=0; %ClimGrid=[-0.05 0.05];   ClimGrid_res=[-0.05 0.05];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch SampleMethod
case {'Quadtree'}
     if SAR_1  charu_1='patch(cx_1,cy_1,datavec_1),                      caxis(ClimQt);title(DataSet_1)';         end
     if SAR_2  charu_2='patch(cx_2,cy_2,datavec_2),                      caxis(ClimQt);title(DataSet_2)';         end
     if SAR_3  charu_3='patch(cx_3,cy_3,datavec_3),                      caxis(ClimQt);title(DataSet_3)';         end
     if SAR_4  charu_4='patch(cx_4,cy_4,datavec_4),                      caxis(ClimQt);title(DataSet_4)';         end
     if SAR_5  charu_5='patch(cx_5,cy_5,datavec_5),                      caxis(ClimQt);title(DataSet_5)';         end

     if SAR_1  charm_1='patch(cx_1,cy_1,transpose(pred(          1:datind(1)))),              caxis(ClimQt)';         end
     if SAR_2  charm_2='patch(cx_2,cy_2,transpose(pred(datind(1)+1:datind(2)))),              caxis(ClimQt)';         end
     if SAR_3  charm_3='patch(cx_3,cy_3,transpose(pred(datind(2)+1:datind(3)))),              caxis(ClimQt)';         end
     if SAR_4  charm_4='patch(cx_4,cy_4,transpose(pred(datind(3)+1:datind(4)))),              caxis(ClimQt)';         end
     if SAR_5  charm_5='patch(cx_5,cy_5,transpose(pred(datind(4)+1:datind(5)))),              caxis(ClimQt)';         end

     if SAR_1  chard_1='patch(cx_1,cy_1,datavec_1-transpose(pred(         1:datind(1)))), caxis(ClimQt_res)';     end
     if SAR_2  chard_2='patch(cx_2,cy_2,datavec_2-transpose(pred(datind(1)+1:datind(2)))),caxis(ClimQt_res)';     end
     if SAR_3  chard_3='patch(cx_3,cy_3,datavec_3-transpose(pred(datind(2)+1:datind(3)))),caxis(ClimQt_res)';     end
     if SAR_4  chard_4='patch(cx_4,cy_4,datavec_4-transpose(pred(datind(3)+1:datind(4)))),caxis(ClimQt_res)';     end
     if SAR_5  chard_5='patch(cx_5,cy_5,datavec_5-transpose(pred(datind(4)+1:datind(5)))),caxis(ClimQt_res)';     end

     if GPShorz charGPShorz   =['H=DM_Quiver(dGPS.xy,dGPS.enu(:),dGPS.dcov,1.0); set(H,''color'',''r''); H=DM_Quiver(mGPS.xy,mGPS.enu(:),mGPS.dcov*0,1.0); set(H,''color'',''b'');']; end
     if GPShorz charGPShorzRes=['H=DM_Quiver(dGPS.xy,dGPS.enu(:)-mGPS.enu(:),dGPS.dcov,1.0); set(H,''color'',''r'')'];  end
     
     if nargout >=1
     if SAR_1 pred_1=griddata(coord_1(1,:)', coord_1(2,:)', pred(          1:datind(1)),xx,yy);                       end
     if SAR_2 pred_2=griddata(coord_2(1,:)', coord_2(2,:)', pred(datind(1)+1:datind(2)),xx,yy);                       end
     if SAR_3 pred_3=griddata(coord_3(1,:)', coord_3(2,:)', pred(datind(2)+1:datind(3)),xx,yy);                       end
     if SAR_4 pred_4=griddata(coord_4(1,:)', coord_4(2,:)', pred(datind(3)+1:datind(4)),xx,yy);                       end
     if SAR_5 pred_5=griddata(coord_5(1,:)', coord_5(2,:)', pred(datind(4)+1:datind(5)),xx,yy);                       end
     end
case {'Grid'}
     if SAR_1 pred_1=griddata(coord_1(1,:)', coord_1(2,:)', pred(          1:datind(1)),xx,yy,'linear');       end  %9/2008: added line,'{'QJ'} to avoid qhull error
     if SAR_2 pred_2=griddata(coord_2(1,:)', coord_2(2,:)', pred(datind(1)+1:datind(2)),xx,yy,'linear');       end  %1/2011 removed {'QJ'} because no long supported in
     if SAR_3 pred_3=griddata(coord_3(1,:)', coord_3(2,:)', pred(datind(2)+1:datind(3)),xx,yy,'linear');       end  % future matlab versions (message in R2010b)
     if SAR_4 pred_4=griddata(coord_4(1,:)', coord_4(2,:)', pred(datind(3)+1:datind(4)),xx,yy,'linear');       end
     if SAR_5 pred_5=griddata(coord_5(1,:)', coord_5(2,:)', pred(datind(4)+1:datind(5)),xx,yy,'linear');       end

     if SAR_1 charu_1='imagesc(x_l,y_l,mod(data_1.*amp_1,modulo)),caxis(ClimGrid),title(DataSet_1)';          end
     if SAR_2 charu_2='imagesc(x_l,y_l,mod(data_2.*amp_2,modulo)),caxis(ClimGrid),title(DataSet_2)';          end
     if SAR_3 charu_3='imagesc(x_l,y_l,mod(data_3.*amp_3,modulo)),caxis(ClimGrid),title(DataSet_3)';          end
     if SAR_4 charu_4='imagesc(x_l,y_l,mod(data_4.*amp_4,modulo)),caxis(ClimGrid),title(DataSet_4)';          end
     if SAR_5 charu_5='imagesc(x_l,y_l,mod(data_5.*amp_5,modulo)),caxis(ClimGrid),title(DataSet_5)';          end

     if SAR_1 charm_1='imagesc(x_l,y_l,mod(pred_1.*amp_1,modulo)),caxis(ClimGrid),title(''Model'')';          end
     if SAR_2 charm_2='imagesc(x_l,y_l,mod(pred_2.*amp_2,modulo)),caxis(ClimGrid),title(''Model'')';          end
     if SAR_3 charm_3='imagesc(x_l,y_l,mod(pred_3.*amp_3,modulo)),caxis(ClimGrid),title(''Model'')';          end
     if SAR_4 charm_4='imagesc(x_l,y_l,mod(pred_4.*amp_4,modulo)),caxis(ClimGrid),title(''Model'')';          end
     if SAR_5 charm_5='imagesc(x_l,y_l,mod(pred_5.*amp_5,modulo)),caxis(ClimGrid),title(''Model'')';          end

     if SAR_1 chard_1='imagesc(x_l,y_l,mod((data_1-pred_1).*amp_1,modulo_res)),caxis(ClimGrid_res)';          end
     if SAR_2 chard_2='imagesc(x_l,y_l,mod((data_2-pred_2).*amp_2,modulo_res)),caxis(ClimGrid_res)';          end
     if SAR_3 chard_3='imagesc(x_l,y_l,mod((data_3-pred_3).*amp_3,modulo_res)),caxis(ClimGrid_res)';          end
     if SAR_4 chard_4='imagesc(x_l,y_l,mod((data_4-pred_4).*amp_4,modulo_res)),caxis(ClimGrid_res)';          end
     if SAR_5 chard_5='imagesc(x_l,y_l,mod((data_5-pred_5).*amp_5,modulo_res)),caxis(ClimGrid_res)';          end

     %%% Yunjun, 2015-12-04: use quiver() instead of DM_Quiver()
     if GPShorz || GPSvert;
         rGPS_enu = dGPS.enu(:)-mGPS.enu(:); rGPS_enu = reshape(rGPS_enu',3,[]);
         %hscale = max(sqrt(sum(dGPS.enu(1:2,:).^2))); vscale = max(dGPS.enu(3,:));
         hscale = 0.5; vscale = 0.5;                                                      % Yunjun, 2015-12-05: temporary for KirishimaPre GPS dataset
         charGPShorz    = ['H=quiver(dGPS.xy(1,:),dGPS.xy(2,:),dGPS.enu(1,:),dGPS.enu(2,:),hscale,''b'',''LineWidth'',2);',...
                           'H=quiver(mGPS.xy(1,:),mGPS.xy(2,:),mGPS.enu(1,:),mGPS.enu(2,:),hscale,''r'');'];
         charGPShorzRes =  'H=quiver(dGPS.xy(1,:),dGPS.xy(2,:),rGPS_enu(1,:),rGPS_enu(2,:),hscale,''r'',''LineWidth'',2);';
         charGPSvert    = ['H=quiver(dGPS.xy(1,:),dGPS.xy(2,:),zeros(size(dGPS.enu(3,:))),dGPS.enu(3,:),vscale,''b'',''LineWidth'',2);',...
                           'H=quiver(mGPS.xy(1,:),mGPS.xy(2,:),zeros(size(mGPS.enu(3,:))),mGPS.enu(3,:),vscale,''r'');'];
         charGPSvertRes =  'H=quiver(dGPS.xy(1,:),dGPS.xy(2,:),zeros(size(rGPS_enu(3,:))),rGPS_enu(3,:),vscale,''r'',''LineWidth'',2);';
         
         % Base autoscale value on average spacing in the x and y
         % directions.  Estimate number of points in each direction as
         % either the size of the input arrays or the effective square
         % spacing if x and y are vectors. (from quiver.m)
         hscale = 1;
         vscale = 1;
         if min(size(dGPS.xy(1,:)))==1, n=sqrt(numel(dGPS.xy(1,:))); m=n; else [m,n]=size(dGPS.xy(1,:)); end
         delx = diff([min(dGPS.xy(1,:)) max(dGPS.xy(1,:))])/n;
         dely = diff([min(dGPS.xy(2,:)) max(dGPS.xy(2,:))])/m;
         del = delx.^2 + dely.^2;
         % Adjust horizontal scale
         if del>0
             len = sqrt((dGPS.enu(1,:).^2 + dGPS.enu(2,:).^2 + dGPS.sig(1,:).^2 + dGPS.sig(2,:).^2)/del);
             maxlen = max(len(:));
         else
             maxlen = 0;
         end
         
         if maxlen>0
             hscale = hscale*0.9 / maxlen;
         else
             hscale = hscale*0.9;
         end
         % Adjust vertial scale
         if del>0
             len = sqrt((zeros(1,length(dGPS.enu(3,:))).^2 + dGPS.enu(3,:).^2 + zeros(1,length(dGPS.enu(3,:))).^2 + dGPS.sig(3,:).^2)/del);
             maxlen = max(len(:));
         else
             maxlen = 0;
         end
         
         if maxlen>0
             vscale = vscale*0.9 / maxlen;
         else
             vscale = vscale*0.9;
         end
         
         %hscale = max((max(east+lons)-min(east+lons))/(x_last-x_first),...
         %             (max(north+lats)-min(north+lats))/(y_first-y_last))*3;
         %hscale = max(sqrt(east.^2+north.^2))*1.2;
         %vscale = max(dGPS.enu(3,:)+dGPS.sig(3,:));
         
     end
%      if GPShorz; charGPS   =['H=DM_Quiver(dGPS.xy,dGPS.enu(:),dGPS.dcov,1.0); set(H,''color'',''b'',''LineWidth'',2); H=DM_Quiver(mGPS.xy,mGPS.enu(:),mGPS.dcov*0,1.0); set(H,''color'',''b'');']; end
%      if GPShorz; charGPSres=['H=DM_Quiver(dGPS.xy,dGPS.enu(:)-mGPS.enu(:),dGPS.dcov,1.0); set(H,''color'',''r'')'];  end

% This is to eliminate weird data points. We should use CLim but this version of PlotModel does not allow this
% Sep 7 2005. Used to be after if Profile. Did not make any sense to be located at that Spot
     if SAR_1; data_1(find(data_1>1))=nan; data_1(find(data_1<-1))=nan; end
     if SAR_2; data_2(find(data_2>1))=nan; data_2(find(data_2<-1))=nan; end   
     if SAR_3; data_3(find(data_3>1))=nan; data_3(find(data_3<-1))=nan; end   
     if SAR_4; data_4(find(data_4>1))=nan; data_4(find(data_4<-1))=nan; end   
     if SAR_5; data_5(find(data_5>1))=nan; data_5(find(data_5<-1))=nan; end   

     if strcmp('2D',Plot)
        G=[radarlook_1;   radarlook_2 ;  radarlook_3  ; radarlook_4   ;radarlook_5   ];
        d=[data_1(:)';data_2(:)';data_3(:)';data_4(:)';data_5(:)'];
        p=[pred_1(:)'    ;pred_2(:)'    ;pred_3(:)'    ;pred_4(:)'    ;pred_5(:)'];
        G(:,2)=0;      % assume u_north=0
        %G(1,:)=[];d(1,:)=[];
        [ux,uy,uz]         =deal(zeros(size(data_1)));
        m=pinv(G)*d;
        [ux(:),uy(:),uz(:)]=deal(m(1,:),m(2,:),m(3,:));
        [data_1,data_2]=deal(ux,uz) ;
        m=pinv(G)*p;
        [ux(:),uy(:),uz(:)]=deal(m(1,:),m(2,:),m(3,:));
        [pred_1,pred_2]=deal(ux,uz) ;
        [SAR_3,SAR_4,SAR_5]=deal(false) ; clear *_d3 *_d4 *_d5
        DataSet_1='east' ; DataSet_2='vertical';
        [vfield(1).data,vfield(2).data]=deal(data_1,data_2) ;
        [vfield(1).data,vfield(2).data]=deal(flipud(data_1),flipud(data_2)) ;
     end
end
if isstruct(Profile) 
   xp           = linspace(Profile.X(1),Profile.X(2),200) ;
   yp           = linspace(Profile.Y(1),Profile.Y(2),200) ;
   profile_dist = sqrt((xp-xp(1)).^2 + (yp-yp(1)).^2);
end

if nargout>=1
     data_mod=[];
     if SAR_1 data_mod(1).data=flipud(data_1);        end
     if SAR_2 data_mod(2).data=flipud(data_2);        end
     if SAR_3 data_mod(3).data=flipud(data_3);        end
     if SAR_4 data_mod(4).data=flipud(data_4);        end
     if SAR_5 data_mod(5).data=flipud(data_5);        end
end
if nargout>=2
     predout=[];
     if SAR_1 predout(1).data=flipud(pred_1.*amp_1); end     % FA 12/2016: included amp_1 to mask
     if SAR_2 predout(2).data=flipud(pred_2.*amp_2); end
     if SAR_3 predout(3).data=flipud(pred_3.*amp_3); end
     if SAR_4 predout(4).data=flipud(pred_4.*amp_4); end
     if SAR_5 predout(5).data=flipud(pred_5.*amp_5); end
end
% TODO: Instead of interpolating u and pred based on dataset(1) coordinates I should recalculate on a dense grid
if nargout>=3
     u=sum(u,2);
    %enu(1).data=griddata(coord_1(1,:)', coord_1(2,:)', u(1:3:3*datind(1)),xx,yy); 
    %enu(2).data=griddata(coord_1(1,:)', coord_1(2,:)', u(2:3:3*datind(1)),xx,yy); 
    %enu(3).data=griddata(coord_1(1,:)', coord_1(2,:)', u(3:3:3*datind(1)),xx,yy); 
    %enu(1).data=flipud(enu(1).data) ;
    %enu(2).data=flipud(enu(2).data) ;
    %enu(3).data=flipud(enu(3).data) ;
     [enu,coord,u]   = ForwardModel_forBasemap(dataset,modelopt,inverseopt.plotdataopt.basemap);
     %[enu,coord,u]   = ForwardModel_forBasemap(modelopt,inverseopt.plotdataopt.basemap,[],[],dataset); %%%%Added ,[],[],dataset %%Anieri 4/27/15
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if Plot==0.5 return ; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigurePaperUnits','centimeters')
set(0,'DefaultFigurePosition',[259 75 768 960],'DefaultFigurePaperPosition',[0.5 1.0 21.0 25.0],...
      'DefaultFigureInvertHardcopy','off','DefaultFigureColor','w','DefaultAxesColor','w','DefaultPatchLineStyle','none');
scsz = get(0,'ScreenSize');

if SAR_1 || SAR_2; fig1=figure('Position',[10 10 0.75*scsz(4) 0.9*scsz(4)],'Name','SAR 1&2 Data and Model'); end;

if SAR_1;

    switch SampleMethod
        case {'Quadtree'}, colormap('default');cm=[ 1 1 1 ; colormap('jet') ; 1 1 1 ] ; colormap(cm);
        case {'Grid'}, colormap('default');q=colormap(jet);cm=[1 1 1; q(33:61,:);q(6:32,:) ] ; colormap(cm);
    end
    
    %%% Plot InSAR Data and GPS Horizontal
    h1=axes('Parent',fig1);
    eval(charu_1); axis xy ; axis image ; grid on;
    hold on; displotmulti(par,objfunc,modelopt,x_unit);
    if GPShorz; eval(charGPShorz); end
    if isstruct(Profile); hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
    if isstruct(Faults);  hold on; plot(Faults.X, Faults.Y, '-w');  hold off; end
    %legend('Source','GPS Data','GPS Model'); legend('boxoff');
    axis equal; axis tight;
    set(h1,'NextPlot','add', 'position',[0.1 0.72 0.32 0.25],'Color',[0 0 0],'XLim',[min(x_l),max(x_l)],'YLim',[min(y_l),max(y_l)])
    
    %%% Plot Model and GPS Vertical
    h1=axes('Parent',fig1);
    eval(charm_1); axis xy; axis image ;    grid on;
    hold on; displotmulti(par,objfunc,modelopt,x_unit);
    if     GPSvert; eval(charGPSvert);
    elseif GPShorz; eval(charGPShorz);
    end
    % if GPShorz; eval(charGPShorz); end
    if isstruct(Profile); hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
    if isstruct(Faults);  hold on; plot(Faults.X, Faults.Y, '-w'); hold off; end
    axis equal; axis tight;
    set(h1,'NextPlot','add','position',[0.1 0.46 0.32 0.25],'xticklabel','','Color',[0 0 0],'XLim',[min(x_l),max(x_l)],'YLim',[min(y_l),max(y_l)])

    %%% Plot Residual or Profile
    h1=axes('Parent',fig1);
    if ~isstruct(Profile)
        eval(chard_1); axis xy; axis image
        hold on; displotmulti(par,objfunc,modelopt,x_unit);   hold off; grid on
        if isstruct(Profile) hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
        if isstruct(Faults)  hold on; plot(Faults.X, Faults.Y, '-w'); hold off; end
        set(h1,'FontSize',12,'NextPlot','add','position',[0.1 0.21 0.32 0.25],'Color',[1 1 1],'XLim',[min(x_l),max(x_l)],'YLim',[min(y_l),max(y_l)]);
    else
        if  strcmp('Grid',SampleMethod)
            d=interp2(xx,yy,data_1,xp,yp,'nearest');
            p=interp2(xx,yy,pred_1,xp,yp,'nearest');
        elseif  strcmp('Quadtree',SampleMethod)
            yc=cy_1(3,:)+(cy_1(3,:)-cy_1(1,:))/2;
            xc=cx_1(1,:)+(cx_1(2,:)-cx_1(1,:))/2;
            d=griddata(xc,yc,datavec_1,xp,yp,'nearest');
            p=griddata(xc,yc,pred(1:datind(1)),xp,yp,'nearest');
        end
        
        if ~ProfileFillWholes; p(find(isnan(d))) = nan; end
        plot(profile_dist,d,'-b',profile_dist,p,'-g','LineWidth',2);
        legend('Data','Model','Location','northeast'); legend('boxoff');
        ylabel(['Motion [' dataset(1).Unit ']']); xlabel(['Distance [' x_unit ']']);
        set(h1,'NextPlot','add','position',[0.1 0.21 0.32 0.25],'Color',[1 1 1]);
    end
    
    % write model
    h1=axes('Parent',fig1);
    text(0,0.5,plotstr,'Units','normalized','FontSize',FontSize,'FontName','Courier')
    set(h1,'NextPlot','add','position',[0.05 0.025 0.32 0.12],'vis','off')
    
    if ischar(out_name); logplot('',[out_name '3'],'') ;  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if SAR_2
    h1=axes('Parent',fig1);
    eval(charu_2); axis xy ; axis image ; grid on
    hold on; displotmulti(par,objfunc,modelopt,x_unit);
    if GPShorz; eval(charGPShorz); end
    if isstruct(Profile) hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
    if isstruct(Faults)  hold on; plot(Faults.X, Faults.Y, '-w'); hold off; end
    set(h1,'NextPlot','add','position',[0.5 0.720 0.32 0.25],'Color',[0 0 0],'XLim',[min(x_l),max(x_l)],'YLim',[min(y_l),max(y_l)]);
    h2=colorbar;
    set(h2,'position',[0.85 0.720 0.03 0.25]);
    
    h1=axes('Parent',fig1);
    eval(charm_2); axis xy; axis image
    hold on; displotmulti(par,objfunc,modelopt,x_unit); grid on
    if     GPSvert; eval(charGPSvert);
    elseif GPShorz; eval(charGPShorz);
    end
    if isstruct(Profile) hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
    if isstruct(Faults)  hold on; plot(Faults.X, Faults.Y, '-w'); hold off; end
    set(h1,'NextPlot','add','position',[0.5 0.46 0.32 0.25],'xticklabel','','yticklabel','','Color',[0 0 0],'XLim',[min(x_l),max(x_l)],'YLim',[min(y_l),max(y_l)])
    h2=colorbar;
    set(h2,'position',[0.85 0.46 0.03 0.25]);
    
    h1=axes('Parent',fig1);
    if ~isstruct(Profile)
        eval(chard_2); axis xy; axis image
        hold on; displotmulti(par,objfunc,modelopt,x_unit);   hold off; grid on
        if isstruct(Profile) hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
        if isstruct(Faults)  hold on; plot(Faults.X, Faults.Y, '-w'); hold off; end
        set(h1,'NextPlot','add','position',[0.5 0.210 0.32 0.25],'yticklabel','','Color','w','XLim',[min(x_l),max(x_l)],'YLim',[min(y_l),max(y_l)])
        h2=colorbar;
        set(h2,'position',[0.85 0.21 0.03 0.25])
    else
        if  strcmp('Grid',SampleMethod)
            d=interp2(xx,yy,data_2,xp,yp,'nearest');
            p=interp2(xx,yy,pred_2,xp,yp,'nearest');
        elseif  strcmp('Quadtree',SampleMethod)
            yc=cy_2(3,:)+(cy_2(3,:)-cy_2(1,:))/2;
            xc=cx_2(1,:)+(cx_2(2,:)-cx_2(1,:))/2;
            d=griddata(xc,yc,datavec_2,xp,yp,'nearest');
            p=griddata(xc,yc,pred(datind(1)+1:datind(2)),xp,yp,'nearest');
        end
        if ~ProfileFillWholes p(find(isnan(d))) = nan; end
        plot(profile_dist,d,'-b',profile_dist,p,'-g','LineWidth',2)  ;
        %ylabel(['motion [' dataset(1).Unit ']']); 
        xlabel(['Distance [' x_unit ']']);
        set(h1,'NextPlot','add','position',[0.5 0.210 0.32 0.25],'Color','w');
    end ;
end
%if (SAR_1 || SAR_2 ) && PostScript print -dpsc a1.ps ;  end

if (SAR_1 || SAR_2 ) && ischar(out_name) logplot('',[out_name '1'],'') ;  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if SAR_3; figure
    
    switch SampleMethod
        case {'Quadtree'},   colormap('default');cm=[ 1 1 1 ; colormap('jet') ; 1 1 1 ] ; colormap(cm);
        case {'Grid'}, colormap('default');q=colormap(jet);cm=[1 1 1; q(33:61,:);q(6:32,:) ] ; colormap(cm);
    end
    
    h1=axes;
    eval(charu_3); axis xy ; axis image ; grid on
    hold on; displotmulti(par,objfunc,modelopt,x_unit);   hold off;
    if isstruct(Profile) hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
    if isstruct(Faults)  hold on; plot(Faults.X, Faults.Y, '-w'); hold off; end
    set(h1,'NextPlot','add', 'position',[0.1 0.700 0.32 0.25],'Color',[0 0 0],'XLim',[min(x_l),max(x_l)],'YLim',[min(y_l),max(y_l)])
    
    axes
    eval(charm_3); axis xy; axis image ;
    hold on; displotmulti(par,objfunc,modelopt,x_unit);   hold off; grid on
    if isstruct(Profile) hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
    if isstruct(Faults)  hold on; plot(Faults.X, Faults.Y, '-w'); hold off; end
    set(gca,'NextPlot','add','position',[0.1 0.43 0.32 0.25],'xticklabel','','Color',[0 0 0],'XLim',[min(x_l),max(x_l)],'YLim',[min(y_l),max(y_l)])
    
    axes
    if ~isstruct(Profile)
        eval(chard_3); axis xy; axis image
        hold on; displotmulti(par,objfunc,modelopt,x_unit);   hold off; grid on
        if isstruct(Profile) hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
        if isstruct(Faults)  hold on; plot(Faults.X, Faults.Y, '-w'); hold off; end
        set(gca,'FontSize',12,'NextPlot','add','position',[0.1 0.160 0.32 0.25],'Color','w','XLim',[min(x_l),max(x_l)],'YLim',[min(y_l),max(y_l)])
    else
        if  strcmp('Grid',SampleMethod)
            d=interp2(xx,yy,data_3,xp,yp,'nearest');
            p=interp2(xx,yy,pred_3,xp,yp,'nearest');
        elseif  strcmp('Quadtree',SampleMethod)
            xc=cx_3(1,:)+(cx_3(2,:)-cx_3(1,:))/2;
            yc=cy_3(3,:)+(cy_3(3,:)-cy_3(1,:))/2;
            d=griddata(xc,yc,datavec_3,xp,yp,'nearest');
            p=griddata(xc,yc,pred(datind(2)+1:datind(3)),xp,yp,'nearest');
        end
        plot(profile_dist,d,'-b',profile_dist,p,'-g','LineWidth',2)  ;
        ylabel(['motion [' dataset(1).Unit ']']); xlabel(['distance[' x_unit ']']);
        set(gca,'NextPlot','add','position',[0.1 0.160 0.32 0.25],'Color','w')
    end ;
    
    % write model
    axes
    text(0,0.5,plotstr,'Units','normalized','FontSize',FontSize,'FontName','Courier')
    set(gca,'NextPlot','add','position',[0.05 0.025 0.32 0.12],'vis','off')
    
end

if SAR_4
    h1=axes;
    eval(charu_4); axis xy ; axis image ; grid on
    hold on; displotmulti(par,objfunc,modelopt,x_unit);   hold off;
    if isstruct(Profile) hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
    if isstruct(Faults)  hold on; plot(Faults.X, Faults.Y, '-w'); hold off; end
    set(h1,'NextPlot','add','position',[0.5 0.700 0.32 0.25],'Color',[0 0 0],'XLim',[min(x_l),max(x_l)],'YLim',[min(y_l),max(y_l)]);
    h2=colorbar;
    set(h2,'position',[0.85 0.700 0.03 0.25]);
    
    axes
    eval(charm_4); axis xy; axis image
    hold on; displotmulti(par,objfunc,modelopt,x_unit);   hold off; grid on
    if isstruct(Profile) hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
    if isstruct(Faults)  hold on; plot(Faults.X, Faults.Y, '-w'); hold off; end
    set(gca,'NextPlot','add','position',[0.5 0.43 0.32 0.25],'xticklabel','','yticklabel','','Color',[0 0 0],'XLim',[min(x_l),max(x_l)],'YLim',[min(y_l),max(y_l)])
    h2=colorbar;
    set(h2,'position',[0.85 0.43 0.03 0.25]);
    
    axes ;
    if ~isstruct(Profile)
        eval(chard_4); axis xy; axis image
        hold on; displotmulti(par,objfunc,modelopt,x_unit);   hold off; grid on
        if isstruct(Profile) hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
        if isstruct(Faults)  hold on; plot(Faults.X, Faults.Y, '-w'); hold off; end
        set(gca,'NextPlot','add','position',[0.5 0.160 0.32 0.25],'yticklabel','','Color','w','XLim',[min(x_l),max(x_l)],'YLim',[min(y_l),max(y_l)])
        h2=colorbar;
        set(h2,'position',[0.85 0.160 0.03 0.25])
    else
        if  strcmp('Grid',SampleMethod)
            d=interp2(xx,yy,data_4,xp,yp,'nearest');
            p=interp2(xx,yy,pred_4,xp,yp,'nearest');
        elseif  strcmp('Quadtree',SampleMethod)
            xc=cx_4(1,:)+(cx_4(2,:)-cx_4(1,:))/2;
            yc=cy_4(3,:)+(cy_4(3,:)-cy_4(1,:))/2;
            d=griddata(xc,yc,datavec_4,xp,yp,'nearest');
            p=griddata(xc,yc,pred(datind(3)+1:datind(4)),xp,yp,'nearest');
        end
        plot(profile_dist,d,'-b',profile_dist,p,'-g')  ;
        ylabel(['motion [' dataset(1).Unit ']']); xlabel(['distance[' x_unit ']']);
        set(gca,'NextPlot','add','position',[0.5 0.160 0.32 0.25],'Color','w');
    end ;
    
end
%if (SAR_3 || SAR_4 ) && PostScript print -dpsc a2.ps ;  end
if (SAR_3 || SAR_4 ) && ischar(out_name) logplot('',[out_name '2'],'') ;  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SAR_5
    figure
    
    switch SampleMethod
        case {'Quadtree'},   colormap('default');cm=[ 1 1 1 ; colormap('jet') ; 1 1 1 ] ; colormap(cm);
        case {'Grid'}, colormap('default');q=colormap(jet);cm=[1 1 1; q(33:61,:);q(6:32,:) ] ; colormap(cm);
    end
    %Plot descending data
    
    h1=axes;
    eval(charu_5); axis xy ; axis image ; grid on
    hold on; displotmulti(par,objfunc,modelopt,x_unit);   hold off;
    if isstruct(Profile) hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
    if isstruct(Faults)  hold on; plot(Faults.X, Faults.Y, '-w'); hold off; end
    set(h1,'NextPlot','add', 'position',[0.1 0.700 0.32 0.25],'Color',[0 0 0])
    
    axes
    eval(charm_5); axis xy; axis image ;
    hold on; displotmulti(par,objfunc,modelopt,x_unit);   hold off; grid on
    if isstruct(Profile) hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
    if isstruct(Faults)  hold on; plot(Faults.X, Faults.Y, '-w'); hold off; end
    set(gca,'NextPlot','add','position',[0.1 0.43 0.32 0.25],'xticklabel','','Color',[0 0 0],'XLim',[min(x_l),max(x_l)],'YLim',[min(y_l),max(y_l)])
    
    axes
    if ~isstruct(Profile)
        eval(chard_5); axis xy; axis image
        hold on; displotmulti(par,objfunc,modelopt,x_unit);   hold off; grid on
        if isstruct(Profile) hold on; plot(Profile.X,Profile.Y,'--k'); hold off; end
        if isstruct(Faults)  hold on; plot(Faults.X, Faults.Y, '-w'); hold off; end
        set(gca,'FontSize',12,'NextPlot','add','position',[0.1 0.160 0.32 0.25],'Color','w','XLim',[min(x_l),max(x_l)],'YLim',[min(y_l),max(y_l)])
    else
        if  strcmp('Grid',SampleMethod)
            d=interp2(xx,yy,data_5,xp,yp,'nearest');
            p=interp2(xx,yy,pred_5,xp,yp,'nearest');
        elseif  strcmp('Quadtree',SampleMethod)
            xc=cx_5(1,:)+(cx_5(2,:)-cx_5(1,:))/2;
            yc=cy_5(3,:)+(cy_5(3,:)-cy_5(1,:))/2;
            d=griddata(xc,yc,datavec_5,xp,yp,'nearest');
            p=griddata(xc,yc,pred(datind(4)+1:datind(5)),xp,yp,'nearest');
        end
        plot(profile_dist,d,'-b',profile_dist,p,'-g')  ;
        ylabel(['motion [' dataset(1).Unit ']']); xlabel(['distance[' x_unit ']']);
        set(gca,'NextPlot','add','position',[0.1 0.160 0.32 0.25],'Color','w')
    end ;
    
    % write model
    % write model
    axes
    text(0,0.5,plotstr,'Units','normalized','FontSize',FontSize,'FontName','Courier')
    set(gca,'NextPlot','add','position',[0.05 0.025 0.32 0.12],'vis','off')
    
    %if PostScript print -dpsc a3.ps ;  end
    if ischar(out_name) logplot('',[out_name '3'],'') ;  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if GPShorz || GPSvert
    
    %%% Yunjun, 2015-12-03: Fix figure size, move subplot2 a little bit up
    %%% to avoid overlap with txt
    scsz = get(0,'ScreenSize');
    figure('Position',[10 10 0.5*scsz(4) 0.9*scsz(4)],'Name','GPS Data and Model');

    %%% To plot DEM shade on the background
    bmapDir = strrep(out_name,'/DataModelResidual_local','');
    load([bmapDir,'/basemap.mat']);
    crnr_lola = [basemap.x_first,basemap.x_first+basemap.x_step*basemap.width;...
                 basemap.y_first,basemap.y_first+basemap.y_step*basemap.file_length];
    crnr_xy = lola2xy(crnr_lola,basemap,1);

    if GPShorz
        if GPShorz && GPSvert; h1=subplot(2,1,1); end;
        %H=DM_Quiver(dGPS.xy,dGPS.enu(:),dGPS.dcov,1.0) ; set(H,'color','r')
        %H=DM_Quiver(mGPS.xy,mGPS.enu(:),mGPS.dcov*0,1.0) ; set(H,'color','b')   % FA 8/08 commented out. Did not work  for Gina's ECSZ data
        image(crnr_xy(1,:),crnr_xy(2,:),basemap.shade);
        hold on ;
        quiver(dGPS.xy(1,:),dGPS.xy(2,:),dGPS.enu(1,:)*hscale,dGPS.enu(2,:)*hscale,0,'b','LineWidth',2)
        ellipse(dGPS.sig(1,:)*hscale,dGPS.sig(2,:)*hscale,zeros(size(dGPS.enu(1,:))),dGPS.xy(1,:)+dGPS.enu(1,:)*hscale,dGPS.xy(2,:)+dGPS.enu(2,:)*hscale,'b')
        hold on ; quiver(mGPS.xy(1,:),mGPS.xy(2,:),mGPS.enu(1,:)*hscale,mGPS.enu(2,:)*hscale,0,'r')
        
        hold on; displotmulti(par,objfunc,modelopt,x_unit);   hold off;
        
        % FA 8/08: need to convert into Faults into xy
        %if plotdataopt.Faults(1)    hold on ; plot(plotdataopt.Faults(:,1),plotdataopt.Faults(:,2),'b') ; hold off ; end
        %if exist('Faults')    hold on ; plot(Faults(:,1),Faults(:,2),'b') ; hold off ; end
        if exist('Sites') && isstruct(Sites); hold on ; text  (Sites.xy(1,:),Sites.xy(2,:),Sites.sites(:,:),'FontSize',6,'VerticalAlignment','bottom','Color','k')   ;hold off;end
        title('Horizontal deformation'); box on;
        legend('Data','Model','Source'); legend('boxoff')
        axis equal; axis tight;
        set(gca,'YDir','normal')
        %         xlim = get(h1,'Xlim'); xwid = xlim(2)-xlim(1);  xlim(xlim+[xwid*-0.05 xwid*0.05]);
        %         ylim = get(h1,'Ylim'); ywid = ylim(2)-ylim(1);  ylim(ylim+[ywid*-0.05 ywid*0.05]);
    end
    if GPSvert
        if GPShorz && GPSvert; h2=subplot(2,1,2); pos1 = get(h1,'Position');  set(h2,'Position',[pos1(1) 0.18 pos1(3:4)]);    end;
        denu=zeros(3,length(dGPS.enu)); [denu(1,:),denu(2,:),denu(3,:)]=deal( zeros(1,length(dGPS.enu)) , dGPS.enu(3,:) , zeros(1,length(dGPS.enu)) ) ;
        menu=zeros(3,length(dGPS.enu)); [menu(1,:),menu(2,:),menu(3,:)]=deal( zeros(1,length(dGPS.enu)) , mGPS.enu(3,:) , zeros(1,length(dGPS.enu)) ) ;
        %            H=DM_Quiver(dGPS.xy,denu(:),dGPS.dcov*0,1.0) ; set(H,'color','r')
        %            H=DM_Quiver(mGPS.xy,menu(:),mGPS.dcov*0,1.0) ; set(H,'color','b')
        image(crnr_xy(1,:),crnr_xy(2,:),basemap.shade);
        hold on ; quiver(dGPS.xy(1,:),dGPS.xy(2,:),zeros(1,length(dGPS.enu(3,:))),dGPS.enu(3,:)*vscale,0,'b','LineWidth',2)
        hold on ; errorbar(dGPS.xy(1,:),dGPS.xy(2,:)+dGPS.enu(3,:)*vscale,dGPS.sig(3,:)*vscale,'bo','MarkerSize',3)
        hold on ; quiver(mGPS.xy(1,:),mGPS.xy(2,:),zeros(1,length(mGPS.enu(3,:))),mGPS.enu(3,:)*vscale,0,'r')
        
        hold on; displotmulti(par,objfunc,modelopt,x_unit);   hold off;
        %if exist('Faults')    hold on ; plot(Faults(:,1),Faults(:,2),'b') ; hold off ; end
        if exist('Sites') && isstruct(Sites); hold on ; text  (Sites.xy(1,:),Sites.xy(2,:),Sites.sites(:,:),'FontSize',6,'VerticalAlignment','bottom','Color','k')   ;hold off;end
        title('Vertical deformation'); box on;
        legend('Data','Model','Source'); legend('boxoff')
        axis equal; axis tight;
        set(gca,'YDir','normal')
    end
    if GPShorz && GPSvert;
    xlim = [min([get(h1,'Xlim'),get(h2,'Xlim')]),max([get(h1,'Xlim'),get(h2,'Xlim')])];
    ylim = [min([get(h1,'Ylim'),get(h2,'Ylim')]),max([get(h1,'Ylim'),get(h2,'Ylim')])];
    xwid = xlim(2)-xlim(1);             xlim = xlim+[xwid*-0.05 xwid*0.05];
    ywid = ylim(2)-ylim(1);             ylim = ylim+[ywid*-0.05 ywid*0.05];
    set(h1,'Xlim',xlim,'Ylim',ylim);    set(h2,'Xlim',xlim,'Ylim',ylim)
    end
    % write model
    axes
    text(0,0.5,plotstr,'Units','normalized','FontSize',FontSize,'FontName','Courier')
    set(gca,'NextPlot','add','position',[0.05 0.025 0.32 0.12],'vis','off')
    
    %     if GPShorz && GPSvert; linkaxes([h1,h2],'xy'); end                      % Synchronize limits of two subplots
    
    %if PostScript print -dpsc a4.ps ;  end
    if ischar(out_name); logplot('',[out_name '4'],'') ;  end ;
    
end
end
