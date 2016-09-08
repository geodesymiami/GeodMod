function [inverseopt] = ModifyBoundsAndParNames(dataset,inverseopt,modelopt)
%   ModifyBounds   - Modify bounds array for non-linear inversion for factors for datasets 
%
%   usage: [bounds] = ModifyBoundsAndParNames(bounds,dataset,inverseopt);
%
%
% INPUT:
%		bounds         - existing bounds  
%       inversetopt    - options structure for objective function
%       dataset        - dataset structure
%       FactorNonLinDelta -if given, variation from one for new bounds    [default 0.2]
%
% OUTPUT:
%       bounds         - modified bounds
%      
% FA Dec 2006 modified so that FactorNonLinDelta is given over inverseopt.FactorNonLinDelta
%
f=fieldnames(inverseopt); for i=1:length(f) eval([char(f{i}) '= inverseopt.(f{i}) ;' ]); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%[SAR_1,SAR_2,SAR_3,SAR_4,SAR_5]=deal(false);
%readfrom_dataset_structure ; 
%iSAR=SAR_1+SAR_2+SAR_3+SAR_4+SAR_5;
N_SAR=sum([dataset.SAR]);
N_GPS=sum([dataset.GPS]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if      strcmp('SAR',FactorNonLin) || strcmp('GPS',FactorNonLin)
   bounds(end+1,1)=1-FactorNonLinDelta;
   bounds(end  ,2)=1+FactorNonLinDelta;
elseif  strcmp('SARmul',FactorNonLin)
   %ifac= N_SAR-1 + (GPShorz || GPSvert);
   ifac= N_SAR-1 + logical(N_GPS);
   ind = [ length(bounds)+1:length(bounds)+ifac ];
   bounds(ind,1)=1-FactorNonLinDelta;
   bounds(ind,2)=1+FactorNonLinDelta;
end

% FA Mar08: significantly simplified the following. Hope I did not break anything
[fac_names,fac_namform,fac_varform,fac_freeind]=deal([]);
if inverseopt.FactorNonLin
   if strcmp('SAR',inverseopt.FactorNonLin)
      fac_names    ={'SARfac'};
      fac_namform  ={'%6s'   };
      fac_varform  ={'%6.2f' };
      fac_freeind  =length(freeind)+length(fixind)+length(linearind)+1;
   elseif strcmp('GPS',inverseopt.FactorNonLin)
      fac_names    ={'GPSfac'};
      fac_namform  ={'%6s'   };
      fac_varform  ={'%6.2f' };
      fac_freeind  =length(freeind)+length(fixind)+length(linearind)+1;
   elseif strcmp('SARmul',inverseopt.FactorNonLin)
      [fac_names,fac_namform,fac_varform,fac_freeind]=deal([]);
      for i=1:ifac                                         % FA Mar08 replaced N_SAR-1 by ifac
          fac_names{i}  =['SAR-' dataset(i).DataSet];
          fac_namform{i}=['%6s'];
          fac_varform{i}=['%6.2f'] ;
          fac_freeind(i)=length(freeind)+length(fixind)+length(linearind)+i;
      end
   end
end

inverseopt.bounds       = bounds;
inverseopt.ParNames     = [modelopt.ParNames     fac_names   ];
inverseopt.ParForm      = [modelopt.ParForm      fac_varform ];
inverseopt.ParNamesForm = [modelopt.ParNamesForm fac_namform ];
inverseopt.freeind      = [inverseopt.freeind;   fac_freeind'];
