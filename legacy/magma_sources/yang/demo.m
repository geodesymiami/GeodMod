% This is demonstration of surface deformations due to a prolate
% shperoid of an arbitrary orientation 
% Yuri Fialko (fialko@gps.caltech.edu), 2/8/00

% LOAD INPUTS
clear all;
eps=1e-10;              % small number
load_par;

[ux,uy,uz]=fcn_yangM(as,x(:),y(:),matrl,tp(:));  % calculate surface displ.
ux=reshape(ux,size(x));
uy=reshape(uy,size(x));
uz=reshape(uz,size(x));
uxsub = ux(1:dhx:ndat,1:dhy:mdat);     % subsample horiz. displ. for output
uysub = uy(1:dhx:ndat,1:dhy:mdat);

%plot model
fig1 = figure(1);
clf
set(fig1,'Name','Surface deformation',...
    'PaperPosition',[0.25 1.5 6 6],...
    'Position',[10 20 600 600],'number','off');
   cmin=min(min(uz));
   cmax=max(max(uz));
   caxis([cmin,cmax])
   axis('equal')
   pcolor(x,y,uz), shading flat, colormap(jet), hold on
   quiver(xsub,ysub,uxsub,uysub,'k'), hold on
   xlabel('W-E distance, km')
   ylabel('S-N distance, km')
   title(['Vert. (color) and Horiz. (arrows) displacements, cm']);
   colorbar('h')
hold off
