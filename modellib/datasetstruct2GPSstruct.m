function [GPSdata,GPSpred]=datasetstruct2GPSstruct(dataset,basemap)
%
% datasetstructure2data   -  extracts GPS data from dataset structure (uses pred as data if given)
%
% usage: function [GPS]=datasetstructure2data(dataset,pred);
%
% Input:   dataset      dataset structure
%
% Output:  GPS structures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ieast   = strmatch('GPSeast' ,{dataset(:).DataSet});
inorth  = strmatch('GPSnorth',{dataset(:).DataSet});
iup     = strmatch('GPSvert' ,{dataset(:).DataSet});

if     ieast;   GPSdata.xy = dataset(ieast).coord;  GPSdata.dcov = dataset(ieast).dcov;
elseif inorth;  GPSdata.xy = dataset(inorth).coord; GPSdata.dcov = dataset(inorth).dcov;
elseif iup;     GPSdata.xy = dataset(iup).coord;    GPSdata.dcov = dataset(iup).dcov;
end

GPSdata.enu = zeros(3,length(GPSdata.xy));
GPSdata.sig = zeros(3,length(GPSdata.xy));

if exist('basemap','var')
    GPSdata.lola = lola2xy(GPSdata.xy,basemap,-1);
end

if ieast;  GPSdata.sig(1,:) = dataset(ieast ).normalization; end
if inorth; GPSdata.sig(2,:) = dataset(inorth).normalization; end
if iup;    GPSdata.sig(3,:) = dataset(iup   ).normalization; end

if ieast;  GPSdata.enu(1,:) = dataset(ieast ).datavec;   end
if inorth; GPSdata.enu(2,:) = dataset(inorth).datavec;   end
if iup;    GPSdata.enu(3,:) = dataset(iup   ).datavec;   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          GPSpred          = GPSdata;

if ieast;  GPSpred.enu(1,:) = dataset(ieast ).predvec;       end
if inorth; GPSpred.enu(2,:) = dataset(inorth).predvec;       end
if iup;    GPSpred.enu(3,:) = dataset(iup   ).predvec;       end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
