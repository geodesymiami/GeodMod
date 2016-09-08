function [locind,locindvect,locind_saved,outind]=LL2ind_igram(igram,InvLoc)

% Transform [Lat Lon] data into Indices for a given igram.data
% [locind,locindvect]=LL2ind(igram,InvLoc)

% N.G. November 2005

if size(InvLoc,2)~=2
    InvLoc=InvLoc' ;
end

Lon0=igram(1).x_first;Lat0=igram(1).y_first;
Incx=igram(1).x_step ;Incy=igram(1).y_step;

locind(:,2)=round((InvLoc(:,2)-Lon0)/Incx+1);
locind(:,1)=round((InvLoc(:,1)-Lat0)/Incy+1);

small1=find(locind(:,1)<0);
small2=find(locind(:,2)<0);
big1=find(locind(:,1)>igram.file_length);
big2=find(locind(:,2)>igram.width);

outind=[small1; small2 ;big1; big2];outind=unique(outind);

locind_saved=locind;
locind(outind,:)=[];

locindvect=YX2L(locind(:,1),locind(:,2),igram.file_length);
