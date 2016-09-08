function [resi,pred,u,rms,npar,Wper,rms_unitsig,mlin,dataset_mod] =GenericObjectiveFunction(invpar,dataset,objfuncopt,sqrflag)
%   GenericObjectiveFunction   - calculates difference between data and model predictions 
%
% INPUT:
%		par	- (*,1) vector of modelparameters  (N*10+4 for N dislocations. Last 4 are for Mogi source)
%		data	- (nx1) vector of observations
%		coord	- (2xn) matrix of coordinates
%		radarlook - (3x1) unit vector in the look direction
%		nu	- poisson's ratio
%		normalization - (1x1)
%		fixind	- indices of fixed model parameters
%		fixpar	- the value of the fixed parameters
%		sqrflag - 1 - the norm of the residual vector is returned
%			  0 - the residual vector is returned (for LEASTSQ).
%      
%       objfuncopt - option for objective function:
%       modelopt         - structure with model specifications
%       PhaseRamp        - linear inversion for PhaseRamp or constant factor ('PhaseRamp,'Const','off')  
%       FactorLin        - linear inversion for strength of sources in individual interferograms             [default 'off']
%       FactorNonLin     - non-linear inversion for factors for individual datasets                          [default 'off']
%                          'SAR'    - one factor for SAR, GPS needs to be given
%                          'SARmul' - multiple factors for SAR (for each dataset if GPS given)
%                          NOTE: objfunc needs to know this to know which parameters are for Sources. The cleaner way to do this
%                          would be a filed objfuncopt.partype=[d d d d d d d d d d m m m m f f f] or similar
%
% OUTPUT
%   	resi	    - either (nx1) residual vector between the input model and the data, or the norm of the vector
%		pred        - 
%		u           - surface displacement
%		rms         - root mean squares
%		npar        - model parameter including results from linear inversion
%		mlin        - all parameters from linear inversion (21 aug 2005: Not sure whether needed anywhere)
%       Wper        - weigths of different datasets in percent
%       rms_unitsig - rms for unit sigmas
%       dataset_mod - modified dataset with the linear phase ramps, constants, and factors removed
%
%  FA, Oct 12 2006   for N_disloc=0 changed to pred_disloc=[]. This reduces warnings and keeps G small. Seems to work (also check ModifyDataSetLin)
%  FA, Dec 8  2006   changed how FactorNonLin is multiplied.
%  FA, Feb 19 2007   now allows for 2 inflation sources (both need to be Mogi sources)
%  FA, June 2007     renamed from ManyDisclocOneInflation. Now uses u=Forwardmodel();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Put in the fixed parameters:
[modelpar]=modelpar2invpar(invpar,objfuncopt,-1);
par       = modelpar ;                                % from here on par is always modelpar (with fixed,linear parameters included)
%RampConstFactorLin = objfuncopt.PhaseRamp(1)+objfuncopt.FactorLin(1) ;   % true if inversion for Ramp,Const or FactorLin   ( prior to feb 08 change)
RampConstFactorLin = objfuncopt.FactorLin(1) + islogical(objfuncopt.PhaseRamp(1))*objfuncopt.PhaseRamp(1) + ischar(objfuncopt.PhaseRamp(1))*strcmp('Const',objfuncopt.PhaseRamp);   % true if inversion for Ramp,Const or FactorLin
if objfuncopt.FactorNonLin [FactorNonLin,sourceend]=GetFactorNonLin(par,dataset,objfuncopt); else [FactorNonLin,sourceend]=deal(false,length(par)); end
[d,coord,normalization,radarlook,datind,hgt,G_phaseramp,D_1,D_2,D_3,D_4,D_5,D_6,D_7,D_8]  =datasetstructure2data(dataset);
[SAR_1,SAR_2,SAR_3,SAR_4,SAR_5,SAR_6,SAR_7,SAR_8]=dealfalk(dataset.SAR); N_SAR=sum([SAR_1,SAR_2,SAR_3,SAR_4,SAR_5,SAR_6,SAR_7,SAR_8]);       % needed for FactorNonLin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% calculate more parameters depending of output parameters for plotting routines
%
    N_disloc         = objfuncopt.modelopt.N_disloc;
    N_mogi           = objfuncopt.modelopt.N_mogi;
    N_penny          = objfuncopt.modelopt.N_penny;
    N_mctigue        = objfuncopt.modelopt.N_mctigue;
    N_yang           = objfuncopt.modelopt.N_yang;
    N_visco1d        = objfuncopt.modelopt.N_visco1d;
    N_lockedandcreep = objfuncopt.modelopt.N_lockedandcreep;
    N_sources        = objfuncopt.modelopt.N_sources;
    linearsources    = objfuncopt.linearsources;
    notlinearsources = objfuncopt.notlinearsources;
    [rngcng_1,rngcng_2,rngcng_3,rngcng_4,rngcng_5,rms,Wper,rms_unitsig]  =deal([]);

    G_pred         = ForwardModel(par,coord,radarlook,objfuncopt.modelopt); % create (Ndata-by-Nsource)-matrix of LOS displacement. G_pred is a Greens Function if the linear parameters are set to 1 and the total
                                                                            % displacement is calculated by G_pred*mlin (mlin determined by linear inversion" mlin=G_pred\d). However, since we allow also for
                                                                            % non-linear inversions for linear parameters (this is historic but unnecessary) we keep those parameters unchanged but set
                                                                            % mlin(nolinearind)=1. We should disallow non-linear inversions for linear parameters (a bit confusing for disloc because 3 parameters can be nonlinear)

if sum(linearsources>0) | RampConstFactorLin
    G              = G_pred;                              % design matrix for all sources
    Glinearsources = G(:,linearsources);                  % design matrix for linear sources only
  
    %
    %  multiply each column of G (i.e. the displacements due to the
    %  individual sources) element-wise with the FactorNonLin vector (Nx1) 

    if FactorNonLin(1) 
       G=G.*repmat(FactorNonLin,1,size(G,2)) ; 
       Glinearsources=Glinearsources.*repmat(FactorNonLin,1,size(G,2)) ; 
    end

   if RampConstFactorLin      
       G=[G G_phaseramp];
       Glinearsources=[Glinearsources G_phaseramp];
   end

%
%  linear inversion for linear model parameters
%
   temp = Glinearsources \ d;                                                                % Does the linear inversion
   
   %Takes model parameters from linear inversion and put them in the right
   %place in mlin.

   if RampConstFactorLin                                                                  % The following was changed Feb 2008 as it did not work for PhaseRamp='Const'
        if islogical(objfuncopt.PhaseRamp(1))*objfuncopt.PhaseRamp(1)                     % for PhaseRamp=true
           mlin=zeros(length(linearsources)+length(notlinearsources)+3*N_SAR,1);            % initialise mlin
           mlin(linearsources,1)    = temp(1:end-3*N_SAR);                                % for linearsources, take results from linear inversion.
           mlin(notlinearsources,1) = 1;                                                  % for parameters that are not linearly inverted, set scale factor to 1 (use input)
           mlin((end-3*N_SAR+1):end)=temp((end-3*N_SAR+1):end);                                        % last lines are ramp paramters - take from linear inversion.
        elseif ischar(objfuncopt.PhaseRamp(1))*strcmp('Const',objfuncopt.PhaseRamp)       % for PhaseRamp='Const'
           mlin=zeros(length(linearsources)+length(notlinearsources)+1*N_SAR,1);            % initialise mlin
           mlin(linearsources,1)    = temp(1:end-1*N_SAR);                                
           mlin(notlinearsources,1) = 1;
           mlin((end-1*N_SAR+1):end)= temp((end-1*N_SAR+1):end);       
        else                                                                              % for FactorLin=true      (FA Feb 08, not sure this works same as for Ramp)
           mlin(linearsources)=temp(1:end-3*N_SAR);                                       % FA Feb 08: This does not ever seem to be used. 
           mlin(notlinearsources,1)=1;                                                    % It probably should be  mlin(linearsources,1)=temp... as above
           mlin=[mlin; temp((end-3*N_SAR+1):end)];
           disp('Attention ! Check comment in GenericObjectiveFunction')
        end
   else
        mlin(linearsources,:)=temp;                                                       % for linearsources, take results from linear inversion.   
        mlin(notlinearsources,:)=1;                                                       % for parameters that are not linearly inverted, set scale factor to 1 (use input)
   end

   pred=G*mlin;                                                                           % makes the forward model for all sources(linear, not linear, ramp)
     
else                                                                                      % makes forward model if there is no need for a linear inversion. 
   mlin=ones(N_sources,1);                                                                % mlin is set to 1's because of multiplication with source strength below

end
     pred=sum(pred,2);                                                                    % sums the LOS changes for all sources
%  calculate the (weighted) difference between data and model

    resi = (pred-d) ./ sqrt(normalization);
    if sqrflag   resi=resi'*resi;   end         % Residual vector or sum or squares?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Return here in calls from inversion routines         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  nargout < 4  return;  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        %

%
% Calculate the regular RMS to be used in plotting program
%
    if  nargout >= 4 
                i=1; ind=[1:datind(1)];              res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res/length(ind)) ;
        if D_2  i=2; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res/length(ind)) ; end
        if D_3  i=3; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res/length(ind)) ; end
        if D_4  i=4; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res/length(ind)) ; end
        if D_5  i=5; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res/length(ind)) ; end
        if D_6  i=6; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res/length(ind)) ; end
        if D_7  i=7; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res/length(ind)) ; end
        if D_8  i=8; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res/length(ind)) ; end
                i=1; ind=[1:datind(1)];              res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res) ;
        if D_2  i=2; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res) ; end
        if D_3  i=3; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res) ; end
        if D_4  i=4; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res) ; end
        if D_5  i=5; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res) ; end
        if D_6  i=6; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res) ; end
        if D_7  i=7; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res) ; end
        if D_8  i=8; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ./sqrt(normalization(ind)); rms(i)=sqrt(res'*res) ; end
% Calculate the model parameters  

if  nargout >= 5
        tmlin=mlin;
        tpar = par;
        npar = [];
        for n=1:N_disloc
            if N_disloc >=n   tpar(8:10) = tpar(8:10)*tmlin(1); npar = [npar; tpar(1:10)]; tpar(1:10)=[]; tmlin(1)=[]; end
        end
        
        if N_mogi   >=1   tpar(4)    = tpar(4)   *tmlin(1); npar = [npar; tpar(1:4) ]; tpar(1:4) =[]; tmlin(1)=[]; end 
        if N_mogi   >=2   tpar(4)    = tpar(4)   *tmlin(1); npar = [npar; tpar(1:4) ]; tpar(1:4) =[]; tmlin(1)=[]; end 
        if N_mogi   >=3   tpar(4)    = tpar(4)   *tmlin(1); npar = [npar; tpar(1:4) ]; tpar(1:4) =[]; tmlin(1)=[]; end 
        if N_penny  >=1   tpar(5)    = tpar(5)   *tmlin(1); npar = [npar; tpar(1:5) ]; tpar(1:5) =[]; tmlin(1)=[]; end 
        if N_mctigue>=1   tpar(5)    = tpar(5)   *tmlin(1); npar = [npar; tpar(1:5) ]; tpar(1:5) =[]; tmlin(1)=[]; end 
        if N_yang   >=1   tpar(4)    = tpar(4)   *tmlin(1); npar = [npar; tpar(1:8) ]; tpar(1:8) =[]; tmlin(1)=[]; end 
        if N_lockedandcreep>=1 tpar(5)=tpar(5)   *tmlin(1); npar = [npar; tpar(1:5) ]; tpar(1:5) =[]; tmlin(1)=[]; end 
        %TODO: InitializeModelopt should generate  modelopt.LinParInd with the index of the linear model parameters.  With this information 
        %      the linear model parameters could be estimated
        %      automatically, without info about the model type (problem: disloc has 3 linear pars)
end
% Calculate the Weights for different data sets
    if  nargout >= 6
        W=1./sqrt(normalization);
                 i=1; ind=[            1:datind(i)];  Wper = [       sum(W(ind))];
        if D_2   i=2; ind=[datind(i-1)+1:datind(i)];  Wper = [Wper ; sum(W(ind))];  end ;
        if D_3   i=3; ind=[datind(i-1)+1:datind(i)];  Wper = [Wper ; sum(W(ind))];  end ;
        if D_4   i=4; ind=[datind(i-1)+1:datind(i)];  Wper = [Wper ; sum(W(ind))];  end ;
        if D_5   i=5; ind=[datind(i-1)+1:datind(i)];  Wper = [Wper ; sum(W(ind))];  end ;
        if D_6   i=6; ind=[datind(i-1)+1:datind(i)];  Wper = [Wper ; sum(W(ind))];  end ;
        if D_7   i=7; ind=[datind(i-1)+1:datind(i)];  Wper = [Wper ; sum(W(ind))];  end ;
        if D_8   i=8; ind=[datind(i-1)+1:datind(i)];  Wper = [Wper ; sum(W(ind))];  end ;
        Wper = Wper/sum(Wper)*100;
    end
% Calculate the RMS for unit sigmas
    if  nargout >= 7 
                i=1; ind=[1:datind(1)];              res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ;
        if D_2  i=2; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ; end
        if D_3  i=3; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ; end
        if D_4  i=4; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ; end
        if D_5  i=5; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ; end
        if D_6  i=6; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ; end
        if D_7  i=7; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ; end
        if D_8  i=8; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ; end
    end
% Return mlin  
    if  nargout >= 8 
                   % nothing needs to be calcualted in this case
    end
% Remove the planes from the dataset structure
    if  nargout >= 9 
               [dataset]=ModifyDatasetLin(dataset,objfuncopt,mlin);  

                i=1; ind=[1:datind(1)];              dataset(i).predvec = pred(ind)';
        if D_2  i=2; ind=[datind(i-1)+1:datind(i)];  dataset(i).predvec = pred(ind)'; end
        if D_3  i=3; ind=[datind(i-1)+1:datind(i)];  dataset(i).predvec = pred(ind)'; end
        if D_4  i=4; ind=[datind(i-1)+1:datind(i)];  dataset(i).predvec = pred(ind)'; end
        if D_5  i=5; ind=[datind(i-1)+1:datind(i)];  dataset(i).predvec = pred(ind)'; end
        if D_6  i=6; ind=[datind(i-1)+1:datind(i)];  dataset(i).predvec = pred(ind)'; end
        if D_7  i=7; ind=[datind(i-1)+1:datind(i)];  dataset(i).predvec = pred(ind)'; end
        if D_8  i=8; ind=[datind(i-1)+1:datind(i)];  dataset(i).predvec = pred(ind)'; end
        dataset_mod=dataset;
    end
    end
