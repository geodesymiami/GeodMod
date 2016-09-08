function []=PlotDistribModel(modelopt,inverseopt)
%PlotDistribModel -  plots east,north,up and up with vectors of model prediction
%
% Falk Amelung, September 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dir_out
opt=modelopt;   f=fieldnames(opt);   for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
if sum(inverseopt.distribopt.viewdir) 
    viewdir = inverseopt.distribopt.viewdir
else
    viewdir = [par.xy(5)+10 20];              % default azimuth clode to strike
end

[s_ss,s_ds,s_up] = deal([]);

% FA 2/2010: reverse negative slip for plotting purposes (needed for Sjonni's Haiti slip distribution) 
% convert par into the patchmodel format if necessary
if size(modelopt.par.xy,2) == 10                             % already in N-by-10 format
       pm       = modelopt.par.xy ;
elseif mod(size(modelopt.par.xy,1),10)==0 && size(modelopt.par.xy,2)==1    % in 10*N-by-1 format
       N_disloc = size(modelopt.par.xy,1)/10;   
       pm       = reshape(modelopt.par.xy',10,N_disloc)';
end
if min(pm(:,8)) < 0 pm(:,8)=-pm(:,8); end
if min(pm(:,9)) < 0 pm(:,9)=-pm(:,9); end
tmppar          = zeros(size(pm,1)*10,1);
tmppar(:)       = [pm]';
modelopt.par.xy = tmppar;
% FA 2/2010: end of sign adjusting code


if sum(distribopt.slip)==1
   if distribopt.slip(1) h(1)=subplot(1,1,1);                          s_ss =plot_disloc3d(modelopt.par.xy,'ss'); end
   if distribopt.slip(2) h(1)=subplot(1,1,sum(distribopt.slip(1:2)));  s_ds =plot_disloc3d(modelopt.par.xy,'ds'); end 
   if distribopt.slip(3) h(1)=subplot(1,1,1);                          s_ds =plot_disloc3d(modelopt.par.xy,'ds'); end 
else
   h(1)=subplot(3,1,1); s_ss  = plot_disloc3d(modelopt.par.xy,'ss');
   h(2)=subplot(3,1,2); s_ds  = plot_disloc3d(modelopt.par.xy,'ds'); 
   h(3)=subplot(3,1,3); s_mag = plot_disloc3d(modelopt.par.xy,'mag');
end

sLim = [0 ceil( max([max(s_ss) max(s_ds) max(s_up)]))];

for i=1:length(h)
    caxis(h(i),sLim);
    xlabel(h(i),'East [km]'); ylabel(h(i),'North [km]'); zlabel(h(i),'Depth [km]');
    view(h(i),viewdir);
end

load slipcol;
cmap = colormap( [ones(round(64*distribopt.PlotThresh/100),3) ; slipcol] );        %PlotThresh is given in percent
colormap(cmap);
slipcolorbaropt=struct(                                   ...
        'Location'       ,        'OutsideLowerRight',    ...
        'colsat'         ,        false             ,     ...
        'width'          ,        0.16              ,     ...
        'height'         ,        0.03              ,     ...
        'x_delta'        ,        0.03              ,     ...
        'y_delta'        ,        0.04              ,     ...
        'label_number'   ,         5                ,     ...
        'Title'          ,        ['slip ''[m]']    ,     ...
        'Cmap'           ,        cmap           ,     ...
        'data'           ,        sLim              ); 
 add_colorbar(slipcolorbaropt);
 
