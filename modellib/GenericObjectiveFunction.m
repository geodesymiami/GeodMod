function [resi,pred,u,rms,npar,Wper,rms_unitsig,mlin,dataset_mod] = ...
    GenericObjectiveFunction(invpar,dataset,objfuncopt,sqrflag)
%% GenericObjectiveFunction - calculates difference between data and model predictions
%
% INPUT:
%		par             - (*,1) vector of modelparameters  (N*10+4 for N dislocations. Last 4 are for Mogi source)
%		data            - (nx1) vector of observations
%		coord           - (2xn) matrix of coordinates
%		radarlook       - (3x1) unit vector in the look direction
%		nu              - poisson's ratio
%		normalization   - (1x1)
%		fixind          - indices of fixed model parameters
%		fixpar          - the value of the fixed parameters
%		sqrflag         - 1: the norm of the residual vector (chi_squared is returned)
%                                 0: the residual vector is returned (for lsqnonlin).
%
%       objfuncopt      - option for objective function:
%       modelopt        - structure with model specifications
%       PhaseRamp       - linear inversion for PhaseRamp or constant factor ('PhaseRamp,'Const','off')
%       FactorLin       - linear inversion for strength of sources in individual interferograms             [default 'off']
%       FactorNonLin    - non-linear inversion for factors for individual datasets                          [default 'off']
%                         'SAR'    - one factor for SAR, GPS needs to be given
%                         'SARmul' - multiple factors for SAR (for each dataset if GPS given)
%                         NOTE: objfunc needs to know this to know which parameters are for Sources. The cleaner way to do this
%                         would be a filed objfuncopt.partype=[d d d d d d d d d d m m m m f f f] or similar
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
%  FA, Jun    2007   renamed from ManyDisclocOneInflation. Now uses u=Forwardmodel();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main Function
%disp('GenericObjective Function called')
[d,coord,sigma,radarlook,datind,hgt,G_phaseramp,...
    D_1,D_2,D_3,D_4,D_5,D_6,D_7,D_8] = datasetstructure2data(dataset);
[G_linearsources,G_nonlinearsources] = deal([]);

[modelpar]      = modelpar2invpar(invpar,objfuncopt,-1);                    % Put in the fixed parameters:
par             = modelpar ;                                                % from here on par is always modelpar (with fixed,linear parameters included)

N_sources       = objfuncopt.modelopt.N_sources;
linearsourceind = objfuncopt.linearsourceind;
N_linearsources = sum(linearsourceind);

% create (Ndata-by-Nsource)-matrix of LOS displacement. G_sources is a Greens Function if the linear parameters are set to 1 and
% the total displacement is calculated by G_sources*m_sources (m_sources=mlin and determined by linear inversion mlin=G_sources\d,
% if all sources are linear). However, since we allow also for non-linear inversions for linear parameters (this is historic but unnecessary)
% we keep those parameters unchanged but set msources(nolinearind)=1. We should disallow non-linear inversions for linear parameters (a bit
% confusing for disloc because 3 parameters can be nonlinear)
[G_sources, u] = ForwardModel(par,coord,radarlook,objfuncopt.modelopt,hgt);

G_linearsources    = [ G_sources(:, linearsourceind)];                      % design matrix for linear sources only (linearsources can be empty)
G_nonlinearsources = [ G_sources(:,~linearsourceind)];                      % design matrix for nonlienar sources

%%%%%%%%%%% Linear Inversion %%%%%%%%%%%%%%
%%% Yunjun, 2015-12-03: Using weight in the same way as dMODELS.
stdRW = 0.0005;                                                             % random walk noise StDev for leveling [0.5 mm/sqrt(yr)], from dMODELS
yrRW = 1;                                                                   % yrw Time interval in years, from dMODELS
SigRW_squared = stdRW^2*yrRW.*ones(size(sigma));
Sig_squared = sigma.^2 + SigRW_squared;
%Sig = sigma.^2 ;                                                           % FA 5/2017: This would be without random walk component
weight = 1./Sig_squared;

    % FA 5/2017: the following lines are from dmodels. I verified that weight=WH.
    % SigWN = diag(sigma.^2);                                                     % diagonal covariance matrix, m^2
    % SigRW = stdRW^2*yrRW*eye(size(SigWN));                                         % random walk noise in m^2, (0.5 mm/sqrt(yr))
    % SigH  = SigWN + SigRW;   
    % WH = SigH\eye(size(SigH));                                                  % weight for least square, 1/m^2 (equivalent to inv(SigH))

mlin = [G_linearsources G_phaseramp] \ ...
       ((d - sum(G_nonlinearsources,2)).*weight);                 %FA 5/2017: keep as 'weight' instead of 'sqrt(weight)' as squaring operation  (e.resi'*resi) will not be performed
% mlin = [G_linearsources G_phaseramp] \ ...
%     ((d - sum(G_nonlinearsources,2))./sqrt(normalization));    %FA 2/2010  takes weighting into account. Important if different datasets have different weight. Not fully verified - need to check equation !
%     mlin = [G_linearsources G_phaseramp] \ ...
%           (d - sum(G_nonlinearsources,2));

m_sources                  = ones(N_sources,1);                             % initiate to 1 for the case of nonlinear sources (e.g.m_sources=[1 1 1]);
m_sources(linearsourceind) = mlin(1:N_linearsources);                       % fill in source parameters from linear inversion (if first of 3 sources is non-linear linearsourceind=[0 1 1]. then e.g m_sources=[1 0.2 -0.3])
m_phaseramp                = mlin(N_linearsources+1:end);
m                          = [m_sources; m_phaseramp];                      % fill in PhaseRamp parameters

if isempty(G_phaseramp); G_phaseramp = zeros(size(G_linearsources)); end    % Newly add, to fix array size difference between G_sources and G_phaseramp, Yunjun, Oct 30 2015
pred = [ G_sources  G_phaseramp ] * m;                                      % calculates data prediction for all sources (linear, not linear, ramp) (G_phaseramp can be empty)
%pred = [ G_sources ] * m;  
pred = sum(pred,2);                                                         % sums the LOS changes for all sources

%%%%%%%% calculate the weighted residual vector and chi_squared  %%%%%%%
resi = (pred-d) .* sqrt(weight);                                            % FA 5/2017: added sqrt as the sigma.^2 (actually tmp=tmp.^2) in datasetstructure2data.m was removed
if sqrflag;   resi=resi'*resi;   end                                        % FA 5/2017: identical to X2 (chi_square) in dModels) (residual vector for lsqnonlin (sqrflag-false)

    % FA 5/2017: the following lines are from dmodels. resi for sqrflag equals X2: X2=(sqrt(weight)'.*(pred-d)')*((pred-d) .* sqrt(weight));
    % r = data - model;                       % residual
    % X2 = r'*Wd*r;                           % Chi Square - full covariance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Return here in calls from inversion routines         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  nargout < 4;  return;  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the weighted RMS to be used in plotting program
if  nargout >= 4
             i=1; ind=[1:datind(i)];              rms(i)=sqrt( sum( (pred(ind)-d(ind)).*(pred(ind)-d(ind)).*weight(ind)/sum(weight(ind)) ) ); 
    if D_2;  i=2; ind=[datind(i-1)+1:datind(i)];  rms(i)=sqrt( sum( (pred(ind)-d(ind)).*(pred(ind)-d(ind)).*weight(ind)/sum(weight(ind)) ) ); end
    if D_3;  i=3; ind=[datind(i-1)+1:datind(i)];  rms(i)=sqrt( sum( (pred(ind)-d(ind)).*(pred(ind)-d(ind)).*weight(ind)/sum(weight(ind)) ) ); end
    if D_4;  i=4; ind=[datind(i-1)+1:datind(i)];  rms(i)=sqrt( sum( (pred(ind)-d(ind)).*(pred(ind)-d(ind)).*weight(ind)/sum(weight(ind)) ) ); end
    if D_5;  i=5; ind=[datind(i-1)+1:datind(i)];  rms(i)=sqrt( sum( (pred(ind)-d(ind)).*(pred(ind)-d(ind)).*weight(ind)/sum(weight(ind)) ) ); end
    if D_6;  i=6; ind=[datind(i-1)+1:datind(i)];  rms(i)=sqrt( sum( (pred(ind)-d(ind)).*(pred(ind)-d(ind)).*weight(ind)/sum(weight(ind)) ) ); end
    if D_7;  i=7; ind=[datind(i-1)+1:datind(i)];  rms(i)=sqrt( sum( (pred(ind)-d(ind)).*(pred(ind)-d(ind)).*weight(ind)/sum(weight(ind)) ) ); end
    if D_8;  i=8; ind=[datind(i-1)+1:datind(i)];  rms(i)=sqrt( sum( (pred(ind)-d(ind)).*(pred(ind)-d(ind)).*weight(ind)/sum(weight(ind)) ) ); end
    rms(i+1)  = sqrt(sum((pred-d).*(pred-d).*weight/sum(weight)));   %weighted RMSD
end

% Calculate the model parameters
if  nargout >= 5
    N_disloc           = objfuncopt.modelopt.N_disloc         ;
    N_fault            = objfuncopt.modelopt.N_fault          ;
    N_mogi             = objfuncopt.modelopt.N_mogi           ;
    N_penny            = objfuncopt.modelopt.N_penny          ;
    N_mctigue          = objfuncopt.modelopt.N_mctigue        ;
    N_yang             = objfuncopt.modelopt.N_yang           ;
    N_multidisloc      = objfuncopt.modelopt.N_multidisloc    ;
    N_visco1d          = objfuncopt.modelopt.N_visco1d        ;
    N_lockedandcreep   = objfuncopt.modelopt.N_lockedandcreep ;
    N_peas             = objfuncopt.modelopt.N_peas           ;
    N_squaredisloc     = objfuncopt.modelopt.N_squaredisloc   ;
    
    if N_multidisloc ml = 10+length(objfuncopt.modelopt.multidislocopt.ind);end
    
    tmlin = m;
    tpar  = par;
    npar  = [];
    
    for n=1:double(N_disloc)
        if N_disloc >=n      tpar(8:10) = tpar(8:10)*tmlin(1); npar = [npar; tpar(1:10)]; tpar(1:10)=[]; tmlin(1)=[]; end
    end
    if N_fault         >=1;   tpar(9)    = tpar(9)   *tmlin(1); npar = [npar; tpar(1:9) ]; tpar(1:9) =[]; tmlin(1)=[]; end
    if N_fault         >=2;   tpar(9)    = tpar(9)   *tmlin(1); npar = [npar; tpar(1:9) ]; tpar(1:9) =[]; tmlin(1)=[]; end
    if N_fault         >=3;   tpar(9)    = tpar(9)   *tmlin(1); npar = [npar; tpar(1:9) ]; tpar(1:9) =[]; tmlin(1)=[]; end
    if N_fault         >=4;   tpar(9)    = tpar(9)   *tmlin(1); npar = [npar; tpar(1:9) ]; tpar(1:9) =[]; tmlin(1)=[]; end
    if N_fault         >=5;   tpar(9)    = tpar(9)   *tmlin(1); npar = [npar; tpar(1:9) ]; tpar(1:9) =[]; tmlin(1)=[]; end
    
    if N_mogi          >=1;   tpar(4)    = tpar(4)   *tmlin(1); npar = [npar; tpar(1:4) ]; tpar(1:4) =[]; tmlin(1)=[]; end
    if N_mogi          >=2;   tpar(4)    = tpar(4)   *tmlin(1); npar = [npar; tpar(1:4) ]; tpar(1:4) =[]; tmlin(1)=[]; end
    if N_mogi          >=3;   tpar(4)    = tpar(4)   *tmlin(1); npar = [npar; tpar(1:4) ]; tpar(1:4) =[]; tmlin(1)=[]; end
    if N_penny         >=1;   tpar(5)    = tpar(5)   *tmlin(1); npar = [npar; tpar(1:5) ]; tpar(1:5) =[]; tmlin(1)=[]; end
    if N_mctigue       >=1;   tpar(5)    = tpar(5)   *tmlin(1); npar = [npar; tpar(1:5) ]; tpar(1:5) =[]; tmlin(1)=[]; end
    if N_yang          >=1;   tpar(4)    = tpar(4)   *tmlin(1); npar = [npar; tpar(1:8) ]; tpar(1:8) =[]; tmlin(1)=[]; end
    if N_multidisloc   >=1;   tpar(8:10) = tpar(8:10)*tmlin(1); npar = [npar; tpar(1:ml)]; tpar(1:ml)=[]; tmlin(1)=[]; end
    if N_squaredisloc  >=1;   tpar(7:9)  = tpar(7:9) *tmlin(1); npar = [npar; tpar(1:9)];  tpar(1:9) =[]; tmlin(1)=[]; end
    if N_lockedandcreep>=1;   tpar(5)    = tpar(5)   *tmlin(1); npar = [npar; tpar(1:5) ]; tpar(1:5) =[]; tmlin(1)=[]; end
    if N_peas          >=1;   npar = tpar ;  end
    %TODO: InitializeModelopt should generate  modelopt.LinParInd with the index of the linear model parameters.  With this information
    %      the linear model parameters could be estimated
    %      automatically, without info about the model type (problem: disloc has 3 linear pars)
end

% Calculate the Weights for different data sets
if  nargout >= 6
    W=weight;
    i=1; ind=[            1:datind(i)];  Wper = [       sum(W(ind))];
    if D_2;   i=2; ind=[datind(i-1)+1:datind(i)];  Wper = [Wper ; sum(W(ind))];  end ;
    if D_3;   i=3; ind=[datind(i-1)+1:datind(i)];  Wper = [Wper ; sum(W(ind))];  end ;
    if D_4;   i=4; ind=[datind(i-1)+1:datind(i)];  Wper = [Wper ; sum(W(ind))];  end ;
    if D_5;   i=5; ind=[datind(i-1)+1:datind(i)];  Wper = [Wper ; sum(W(ind))];  end ;
    if D_6;   i=6; ind=[datind(i-1)+1:datind(i)];  Wper = [Wper ; sum(W(ind))];  end ;
    if D_7;   i=7; ind=[datind(i-1)+1:datind(i)];  Wper = [Wper ; sum(W(ind))];  end ;
    if D_8;   i=8; ind=[datind(i-1)+1:datind(i)];  Wper = [Wper ; sum(W(ind))];  end ;
    Wper = Wper/sum(Wper)*100;
end

% Calculate the RMS for unit sigmas
if  nargout >= 7
    i=1; ind=[1:datind(1)];                       res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ;
    if D_2;  i=2; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ; end
    if D_3;  i=3; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ; end
    if D_4;  i=4; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ; end
    if D_5;  i=5; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ; end
    if D_6;  i=6; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ; end
    if D_7;  i=7; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ; end
    if D_8;  i=8; ind=[datind(i-1)+1:datind(i)];  res = (pred(ind)-d(ind)) ; rms_unitsig(i)=sqrt(res'*res/length(ind)) ; end
end

% Return mlin
if  nargout >= 8
    % nothing needs to be calcualted in this case
end

% Remove the planes from the dataset structure
if  nargout >= 9
    [dataset]=ModifyDatasetLin(dataset,objfuncopt,m_phaseramp);
    
    i=1; ind=[1:datind(1)];                       dataset(i).predvec = pred(ind)';
    if D_2;  i=2; ind=[datind(i-1)+1:datind(i)];  dataset(i).predvec = pred(ind)'; end
    if D_3;  i=3; ind=[datind(i-1)+1:datind(i)];  dataset(i).predvec = pred(ind)'; end
    if D_4;  i=4; ind=[datind(i-1)+1:datind(i)];  dataset(i).predvec = pred(ind)'; end
    if D_5;  i=5; ind=[datind(i-1)+1:datind(i)];  dataset(i).predvec = pred(ind)'; end
    if D_6;  i=6; ind=[datind(i-1)+1:datind(i)];  dataset(i).predvec = pred(ind)'; end
    if D_7;  i=7; ind=[datind(i-1)+1:datind(i)];  dataset(i).predvec = pred(ind)'; end
    if D_8;  i=8; ind=[datind(i-1)+1:datind(i)];  dataset(i).predvec = pred(ind)'; end
    dataset_mod=dataset;
end
