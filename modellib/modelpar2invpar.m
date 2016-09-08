function [npar] = modelpar2invpar(par,objfuncopt,direction)
%   modelpar2invpar - converts model parameters to inversion parameters amd vice-versa
%
% INPUT: 
%       par  	         - model parameter vector (Nx1) or bounds array (Nx2) 
%       objfunctopt      - option structure for objective function
%       dir              - direction (1: model-->inv; -1: inv->model)          (default: 1)
%
% OUTPUT:
%       npar             -inversion parameter vector with fixed parameters removed (filled) (to call disloc1, mogi, etc) 
%
if ~exist('direction','var') direction =1; end  

pardim=size(par,2);                                  

if pardim==0 pardim=1; end                           % FA 07/08:  makes it working if par is empty 

if direction==1
   fixind        = objfuncopt.fixind ;
   par(fixind) = [];
   npar          = par;
else                                                  % direction=-1
   %  Put in the fixed parameters:
   fixind =objfuncopt.fixind; 
   freeind=objfuncopt.freeind;
   fixpar =objfuncopt.fixpar;
   linearind=objfuncopt.linearind;
   npar=zeros(length(freeind)+length(fixind)+length(linearind),pardim);    %
   %  Put in the varying parameters:
   %
    npar(fixind,:) =repmat(fixpar,1,pardim); 
                                                            % FA Mar08:  Note to Juliet: We may want have one if to test whether we are one or multi-dimensional for speed
   if length(freeind)==size(par,1)                          % FA Mar08:  wouldn't this be clearer as:  if linearind
       npar(freeind,:)  =par; 
       npar(linearind,:)=ones(length(linearind),pardim); 
   else                                                      % FA Jul08: true only if par contains the linear parameters (size(par,1)>length(freeind) from GenerateSummary)
       modeledparind=sort([freeind;linearind]); 
       npar(modeledparind,:)=par;                            % FA Mar08: put in by Juliet.  Is this OK for pardim>1 ??
   end        
end
