function [InSARpoints]=GetInSARvalues(InSAR,points,opt)

%
% Get values of the InSAR or SBAS matrice at given locations (geographic, indices, ...) 
%
% [points]=GetInSARvalues(InSAR,points,opt)
%
% 'InSAR'  : SBAS or InSAR matrice
%
% 'points': structure with: 
%		points.data: 2 by ? matrice with datapoints
%		points.Proj: 'GEO', 'IND', ... depending on 
%
% opt can contain 
%
% 'method'       : how to look for the points:
%		'exact'    : get the value at the exact pixel location (default)
%		'nearest'  : get the value at the nearest pixel	
%		'average'  : average by input number of pixel
%
% 'methodopt'    : window size to use for averaging (structure)
%		methodopt.size 
%		methodopt.unit  (Only works with input size in pixels)
%		methodopt.maxdistance (distance between the location points and the closest InSAR pixel (in kilometers))
%
%
% N. Gourmelen, Feb 2007
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load SBAS-GPS.mat
%for i=1:length(GPS)
%	GPSlat(i,:)=GPS(i).lat(1);GPSlon(i,:)=GPS(i).lon(1);
%end

%[locind,locindvect]=LL2ind_igram(SBASttt(1),[GPSlat GPSlon]);

%opt=struct('InSAR',SBAS,'points',struct('data',locind,'Proj','IND'),'method','nearest');

defaultopt=struct(                                                          ...
 		'method'                        	,        'off'      ,           ...	
        'methodopt'                         ,        'off'      )            ;             
    
    
if ~exist('opt','var')  [opt] = read_options_fromfile([mfilename '.min'],[]); end ;
[opt] = process_defaultoptions(opt,defaultopt);  display(opt)
f     = fieldnames(opt) ;
for i=1:length(f)  eval([char(f{i}) '= opt.(f{i}) ;' ]);  end

if ~method  method = 'exact';  end

%%%%%%%%%%%%%%%%%%%%
% Check input data %
%%%%%%%%%%%%%%%%%%%%

%if exist('SBAS')==2
%	load SBAS;
%end

if isstruct(methodopt)
	if ~isfield(methodopt,'maxdistance')
		methodopt.maxdistance = 5 ;
    end
else
    methodopt.maxdistance = 5 ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if input points are indices %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isstruct(InSAR)~=1)  tmpmat = InSAR ;  else  tmpmat = InSAR(1) ;  end

if (points.Proj=='IND')
	Indices = points.data ;
elseif (points.Proj=='GEO')
	Indices = LL2ind_igram(tmpmat(1),points.data)  ;
else
	error (sprintf('No good projection, exiting')) ;
end

clear tmpmat ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One degree equivalent distance at scene location  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
latcenter = (InSAR(1).y_first(1) + InSAR(1).y_step * InSAR(1).file_length/2)                       ;  
loncenter = (InSAR(1).x_first(1) + InSAR(1).x_step * InSAR(1).width/2)                             ;
latlength = distance(latcenter - 0.5,loncenter,latcenter + 0.5,loncenter,almanac('earth','grs80')) ;
lonlength = distance(latcenter,loncenter - 0.5,latcenter,loncenter + 0.5,almanac('earth','grs80')) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the InSAR value at given location %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(method,'exact'))

	if (isstruct(InSAR)~=1)
		
		InSARpoints.data        = val(InSAR,Indices)    ;
		InSARpoints.method      = 'exact'               ;
		InSARpoints.Indices     = Indices               ;
		InSARpoints.Coordinates = points.data           ;
        
    elseif (length(InSAR)==1)
		
		InSARpoints.data        = val(InSAR.data,Indices)    ;
		InSARpoints.method      = 'exact'               ;
		InSARpoints.Indices     = Indices               ;
		InSARpoints.Coordinates = points.data           ;
	elseif (isstruct(InSAR)==1)

		for ni=1:length(InSAR)
			InSARpoints(ni).data        = val(InSAR(ni).data,Indices(:,1),Indices(:,2)) ;
			InSARpoints(ni).method      = 'exact'                                       ;
			InSARpoints(ni).Indices     = Indices                                       ;
        	InSARpoints(ni).Coordinates = points.data                                   ;
			InSARpoints(ni).dateYears  = InSAR(ni).dateYears                          ;
        end
	end

elseif (strcmp(method,'nearest'))

    [GoodIndL,GoodIndC]=find(isfinite(InSAR(2).data)==1);

    for nj=1:length(InSAR)

        InSARpoints(nj).dateYears = InSAR(nj).dateYears ;
        InSARpoints(nj).method     = 'nearest'            ;

        for ni=1:size(Indices,1)

            MinDistance        = sqrt(((GoodIndL-Indices(ni,1))*InSAR(1).y_step*latlength).^2+((GoodIndC-Indices(ni,2))*InSAR(1).x_step*lonlength).^2);
            MinDistanceIndices = find(MinDistance==min(MinDistance));

            if (min(MinDistance) > methodopt.maxdistance)

                InSARpoints(nj).data(ni)       = NaN;
                InSARpoints(nj).Indices(ni,:)  = Indices(ni,:);

            else

                if (max(size(MinDistanceIndices))>1)

                    InSARpoints(nj).data(ni)      = mean(val(InSAR(nj).data,GoodIndL(MinDistanceIndices),GoodIndC(MinDistanceIndices)));
                    InSARpoints(nj).Indices(ni,:) = [GoodIndL(MinDistanceIndices(1)) GoodIndC(MinDistanceIndices(1))];

                else

                    InSARpoints(nj).data(ni)      = InSAR(nj).data(GoodIndL(MinDistanceIndices),GoodIndC(MinDistanceIndices));
                    InSARpoints(nj).Indices(ni,:) = [GoodIndL(MinDistanceIndices(1)) GoodIndC(MinDistanceIndices(1))];
                end
            end
        end
    end

elseif (strcmp(method,'average'))

    if (isstruct(InSAR)==1)
        
        %% Get pixels within the circular distance
        
        [XX,YY]  = meshgrid ([1:InSAR(1).width],[[1:InSAR(1).file_length]]) ;
        mask_gps = repmat(struct('distance',ones(size(XX))),size(Indices,1),1)  ;
        
        for ni=1:size(Indices,1)
            XX_seed  = XX - Indices(ni,2) ;
            YY_seed  = YY - Indices(ni,1) ;
            mask_gps(ni).distance = sqrt(XX_seed.^2 + YY_seed.^2) ;
        end

        for nj=1:length(InSAR)

            InSARpoints(nj).dateYears = InSAR(nj).dateYears;

            for ni=1:size(Indices,1)
                
                InSARpoints(nj).Indices(ni,:)  = find(mask_gps(ni).distance <= methodopt.maxdistance) ;
                InSARpoints(nj).dataList(ni,:) = InSAR(nj).data(InSARpoints(nj).Indices(ni,:))        ;

                if    sum(isfinite(InSARpoints(nj).dataList(ni,:)))==0  InSARpoints(nj).data(ni)      = NaN;
                else  InSARpoints(nj).data(ni)      = mean(InSARpoints(nj).dataList(ni,isfinite(InSARpoints(nj).dataList(ni,:))));
                end
            end
        end
    end
end














