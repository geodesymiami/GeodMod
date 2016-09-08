function [dataset] = extract_hgt_from_dem(dem,dataset,x_posting,y_posting)
%   extract_hgt_from_dem - assigns dem height to dataset structure 
%
%    usage:  [dataset]=extract_hgt_from_dem(dem,dataset)
%
% INPUT:
%		dem                  - DEM in same co-ordinate system as dataset structure
%       dataset              - data structure
%
% OUTPUT:
%       dataset              - modified dataset structure with field dataset.hgt
%
%

% first detect and remove sites out of the area

% TODO: The following is a bad hack. It puts the hgt to zero if only GPS data are given (x_posting not given) 
%       Need to rewrite this function
if ~isfield(dataset,'x_posting') || isempty(dataset(1).x_posting)
    if ~exist('x_posting','var') || ~exist('y_posting','var')
        for i=1:length(dataset)
            dataset(i).hgt=dataset(i).datavec*0;
        end
        warning('hack in extrachgt_from_dem: GPS height simply put to zero')
        return
    end
end

if ~exist('x_posting','var') || ~exist('y_posting','var')
    x_posting=dataset(1).x_posting;
    y_posting=dataset(1).y_posting;
end

for i=1:length(dataset)
    if ~strcmp(dataset(i).CoordSystem,'ProjectProfile')
        [x,y]=deal(dataset(i).coord(1,:),dataset(i).coord(2,:)) ;
        x_ind=ceil(x/x_posting);
        y_ind=ceil(y/y_posting);
        
        % assign points out of DEM area hgt of nearest pt. (may happen for GPS)
        x_ind(find(x_ind<=0))=1 ; x_ind(find(x_ind>length(dem)))=length(dem);
        y_ind(find(y_ind<=0))=1 ; y_ind(find(y_ind>length(dem)))=length(dem);
        
        for j=1:length(x_ind)
            hgt(j)=  dem(y_ind(j),x_ind(j));
        end
        dataset(i).hgt=hgt;
        clear x_ind y_ind hgt
    end
    
    if isfield(dataset(1),'fulldata')
        fulldata=dataset(1).fulldata  ;
        tmp=cell(length(dataset),1); [tmp{:}]=deal(dataset.hgt); hgt=cell2mat(tmp');
        fulldata.hgt=hgt;
        dataset(1).fulldata=fulldata;
    end
end

