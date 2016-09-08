function [data,GPS]=datasetstructure2datastructure(dataset,basemap);
%
% datasetstructure2data   -  converts dataset structure into data format for modelling       
%
% usage: function [d,coord,normalization,radarlook,datind,hgt,D_1,D_2,D_3,D_4,D_5]=datasetstructure2data(dataset);
%
% Input:   dataset           dislocation geometry
%
% Output:  d  ,coord
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[D_1,D_2,D_3,D_4,D_5,D_6,D_7,D_8]                                 = dealfalk(dataset.exist);
[coord_1,coord_2,coord_3,coord_4,coord_5,coord_6,coord_7,coord_8] = dealfalk(dataset.coord);
[Ndata_1,Ndata_2,Ndata_3,Ndata_4,Ndata_5,Ndata_6,Ndata_7,Ndata_8] = dealfalk(dataset.Ndata);
[SAR_1,SAR_2,SAR_3,SAR_4,SAR_5,SAR_6,SAR_7,SAR_8]                 = dealfalk(dataset.SAR);
[GPS_1,GPS_2,GPS_3,GPS_4,GPS_5,GPS_6,GPS_7,GPS_8]                 = dealfalk(dataset.GPS);
[ind_1,ind_2,ind_3,ind_4,ind_5,ind_6,ind_7,ind_8]                 = dealfalk(dataset.ind);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_SAR = sum([dataset(:).SAR])
N_GPS = sum([dataset(:).GPS])

for i=1:N_SAR
    data(i)= dataset(i) ;
end

GPS.enu(1,:) = dataset(N_SAR+1).datavec(:);
GPS.enu(2,:) = dataset(N_SAR+2).datavec(:);
GPS.xy       = dataset(N_SAR+1).coord;
GPS.lola     = lola2xy(GPS.xy,basemap,-1);
GPS.Unit     = dataset(N_SAR+1).Unit;
if N_GPS==3 GPS.enu(3,:) = dataset(N_SAR+3).datavec(:); end       % if vertical given
