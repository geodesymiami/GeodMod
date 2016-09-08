function [modelopt]=InverseLinDistrib(dataset,inverseopt,modelopt);
%function [dislocations,parplane,parinfls]=DistOp(dg1,dataset,distopopt,pafrinfl);
% InverseLinDistrib-  linear inversion for fault slip and opening 
%

%
% fields of disopopt:                                        
% PatchSize        [length width] of each  patch                    [default [1.0 1.0]]                                        
% ExtendSize       [length width] of extended fault                 [default 'off']                                        
% ExtendTorSurface flag to extend input disloc to surface           [default 'off']                                         
% InverseSign      flag to inverse sign of data (because of fnnls)  [default 'on']                                         
% PhaseRamp        remove ramp or constant ('on,'Const','off')      [default 'on']                                         
% kappa            smoothing parameter                              [default 0.1]                                         
% factor           flag to invert for multiplier                    [default 'off']                                         
% inflation        inflation source  ('Mogi' or 'Penny')            [default 'off']                                         
% plot             plot results                                     [default 'off']                                         
%
% V1.0  Falk Amelung, September 2008   based on Sjonni Jonnson's model_step1,2,3 programs
%
%    weighting definition- - identities with Sjonni's code:
%                            W       == 1./sqrt(normalization)
%                            1./W2   == normalization  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% June 2007: The new code should also use ModifyDatasetLin for the PhaseRamp removal.
%%% This will take care of whether GGPS data exist or not and also of
%%% the FactorLin problem if ever solved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dir_out
defaultopt=struct(                                         ...
        'fname'             ,   'distribmodel'       ,     ...
        'objfunc'           ,   'GenericObjectiveFunction', ...
        'FactorLin'         ,   'off'                ,     ...
        'FactorNonLin'      ,   'off'                ,     ...
        'FactorNonLinDelta' ,    0.2                 ,     ...
        'PhaseRamp'         ,    'on'                ,     ...
        'Inflation'         ,   'off'                ,     ...
        'Factor'            ,   'off'                ,     ...
        'InverseSign'       ,   'off'                ,     ...
        'Plot'              ,   'off'      )         ;
defaultopt.distribopt=struct(                                   ...
        'DoIt'              ,   'off'                ,     ...
        'ExtendToSurface'   ,   'on'                 ,     ...
        'PatchSize'         ,  [ 1.0  1.0]           ,     ...
        'ExtendSize'        ,   'off'                ,     ...
        'slip'              ,  [ 0 0 1 ]             ,     ...
        'kappa'             ,    0.1                 ,     ...
        'viewdir'           ,   'off'                ,     ...
        'algorithm'         ,   'fnnls'              ,     ...
        'LineStyle'         ,   '-'                  ,     ...
        'CLim'              ,   'off'                ,     ...
        'PlotThresh'        ,    5                   ,     ...
        'modelSjonni'        ,   'off'                ,     ...
        'InverseSign'       ,   'off'                )     ;

global dir_out
[inverseopt]=process_defaultoptions(inverseopt,defaultopt);  
f=fieldnames(inverseopt) ; for i=1:length(f) eval([char(f{i}) '= inverseopt.(f{i}) ;' ]) ; end

out_name  = fullfile(dir_out,fname); 
out_name2 = fullfile(dir_out,'model');
inverseopt.plotdataopt.Surface3D = false;   %3D plotting not here because PrepareDataModel_* may be called

if  ~distribopt.DoIt  return; end
if  CheckInOut('',out_name)  load(out_name); return; end     % return if distrib.mat exist

if ~isfield(modelopt.par,'xy') 
   logmessage('no dislocation, no inversion for slip/opening distribution possible')
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize similar to InverseNonLin
  [inverseopt]           = InitializeInverseopt(dataset,inverseopt,modelopt);                               % generates objfuncopt
  [dataset]              = InitializeDataset(dataset,inverseopt) ;                                         %adds G_phaseramp and normalization to dataset
  
  [t_inverseopt]         = inverseopt; t_inverseopt.objfuncopt.PhaseRamp = 1;          % modify inverseopt so that phaseramp_dummy is removed even for PhaseRamp=Const and 'off'
  [t_inverseopt]         = inverseopt;           % FA 2/2010: Don't know why PhaseRamp is set to 1 for all cases. Try without.
  m_phaseramp_dummy      = phaseramp_dummy(dataset,t_inverseopt.objfuncopt);
  [tmpdataset]           = ModifyDatasetLin(dataset,t_inverseopt.objfuncopt,m_phaseramp_dummy,'replace');     % 1/11 FA: I would be surprised if this works. wrong arguments??
  [tmpdataset]           = InitializeDataset(tmpdataset,inverseopt) ;                                         %adds G_phaseramp and normalization to dataset
  
  [d,coord,normalization,radarlook,datind,hgt,G_phaseramp,D_1,D_2,D_3,D_4,D_5,D_6,D_7,D_8]  = datasetstructure2data(tmpdataset,inverseopt.objfuncopt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sum(modelopt.par.xy(8:10:end)<0)  errordlg('strike-slip negative, fnnls not possible, change fault geometry, will use svd'); distribopt.algorithm='svd'; end
if sum(modelopt.par.xy(9:10:end)<0)  errordlg(   'dip-slip negative, fnnls not possible, change fault geometry, will use svd'); distribopt.algorithm='svd'; end
if sum(modelopt.par.xy(10:10:end)<0) errordlg(    'opening negative, fnnls not possible, change fault geometry, will use svd'); distribopt.algorithm='svd'; end

slip              = distribopt.slip
nu                = modelopt.dislocopt.nu;      % Poisson's ratio

%if distribopt.InverseSign  d=-1*d;  end       % changing sign of data

 if modelopt.N_multidisloc
     modelopt.N_disloc = modelopt.multidislocopt.N_disloc;
     modelopt.par.xy   = multidislocpar2dislocpar(modelopt.par.xy,modelopt.multidislocopt,'km');
  end

%% Start: loop over dislocations  
G_prime  = [];
G_smooth = [];
pm_prime = [] ;
for i=1:modelopt.N_disloc

   % Setup fault geometry (build patches from iput dislocations)

   ind                         = (i-1)*10+1:(i-1)*10+10;  
   [pm,nhorz,nvert,ExtendSize] = MakeFaultPatches(modelopt.par.xy(ind),distribopt);

   pm_prime = [pm_prime; pm];      % save patchmodel for all dislocations to build the complete dislocation solution below

   % Construct kernel 
  
     [K1,K2,K3]   = MakeRadarKern(pm, slip, coord, nu, nhorz, nvert, radarlook', length(coord));    %length(coord) is from Sjonni, probably not needed
     K_unweighted = [K1 K2 K3];                       

          % try using ForwardModel: It works fine 
            if sum(slip) == 1                  % case slip = [1 0 0]  or [0 1 0] or [0 0 1]
               modelpar     = reshape([pm repmat(slip,size(pm,1),1)]',size(pm,1)*10,1);
            elseif slip == [1 1 0]
               modelpar     = [reshape([pm repmat([1 0 0],size(pm,1),1)]',size(pm,1)*10,1); reshape([pm repmat([0 1 0],size(pm,1),1)]',size(pm,1)*10,1)];
            end

            [tmpinverseopt,distribmodelopt]   = Solution2ModifiedOptions(size(modelpar,1)/10,modelpar,inverseopt,dataset);
            [G_pred,u]                        = ForwardModel(modelpar,coord,radarlook,distribmodelopt);
            %K_unweighted = G_pred;        % FA 2/2010  ForwardModel produces the same result as MakeRadarKern as it should 
 
     K_weighted   = K_unweighted .* repmat(1./sqrt(normalization),1,size(K_unweighted,2));
  
     G_prime       = [G_prime K_weighted];

     % Construct smoothing matrix (second order finite difference matrix (laplacian)) (G' in eq.4 from Sjonni's paper)

     [Lap]    = MakeSmoothingMatrix(slip, pm, nhorz, nvert, ExtendSize);
                if sum(slip) == 2 Lap = [Lap Lap]; end
     G_smooth = [G_smooth  Lap * distribopt.kappa];

end
%%%%End: loop over dislocations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%
% weight data 
%
  dprime = d./sqrt(normalization);                          % d' in eq.3 
%
% Construct matrizes of inversion problem (eq. 4) (fill data vector with zeros)   
%

  G      = [G_prime;  G_smooth                                                    ];   % right-hand side in eq. 4
  G      = [G         [G_phaseramp ; zeros(size(G_smooth,1),size(G_phaseramp,2))] ];   % add phaseramp parameters to design matrix 
  din    = [dprime;    zeros(size(G_smooth,1),1)                                  ];   % left -hand side in eq. 4 

%
%% solve using fast non-negative least squares
%
  
  if strcmp(distribopt.algorithm,'fnnls')
     logmessage('inversion using fast non-negative least squares...' );
     tic
     sal        = fnnls(G'*G,G'*din);
     t          = toc;
     logmessage(sprintf('Elapsed time for fnnls: %.3f sec',t));
  else strcmp(distribopt.algorithm,'svd')
      logmessage('inversion using singular value decomposition...' );
      sal       = pinv(G)*din;
  end
  s          = sal(1:end-size(G_phaseramp,2));% the slip of all faults
  m_phaseramp= sal(end-size(G_phaseramp,2)+1:end)';
  m_phaseramp
  
  % Calculate the solution roughness
  nel        = nhorz * nvert;
  tmp1       = G_smooth*s;
  roughness  = sum(abs(tmp1))/(2*(nel));
  str        = sprintf('Kappa: %1.3f Mean roughness: %1.2f cm/km',distribopt.kappa,roughness*100); disp(str);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inverseopt.plotdataopt.modelopt  = modelopt ;
inverseopt.objfuncopt.modelopt   = modelopt ;

N_disloc                         = size(pm_prime,1);

% put all patches into one-dimensional vector and prepare summary

modelpar                                        = zeros(size(pm_prime,1)*10,1) ;
u_disloc                                        = zeros(size(pm_prime,1),3);
u_disloc(find(repmat(slip,size(pm_prime,1),1))) = s;
modelpar(:)                                     = [pm_prime u_disloc]' ;

if isfield(inverseopt.distribopt,'modelSjonni') && inverseopt.distribopt.modelSjonni(1)
   logmessage('####Using Sjonnis solution####')
   load(inverseopt.distribopt.modelSjonni);
   modtot(:,7)                                  = size(dataset(1).data,1)*dataset(1).y_posting+modtot(:,7);
   N_disloc                                     = size(modtot,1);
   modelpar                                     = zeros(size(modtot,1)*10,1) ;
   modelpar(:)                                  = [modtot]';
end

[tmpinverseopt,distribmodelopt]                 = Solution2ModifiedOptions(N_disloc,modelpar,inverseopt,dataset,roughness);
disp( GenerateSummary(distribmodelopt,dataset,tmpinverseopt) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  Plot and save  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distribmodelopt.distribopt       = distribopt;
modelopt.distribmodelopt         = distribmodelopt;
save(out_name, 'distribmodelopt')
save(out_name2,'modelopt')

tmpinverseopt.distribopt.viewdir = [90-distribmodelopt.par.xy(5) 0];
PlotDistribModel(distribmodelopt,tmpinverseopt);

tmpinverseopt.distribopt.viewdir = false ;
PlotDistribModel(distribmodelopt,tmpinverseopt);

PrepareDataModelResidual_local_AndPlot(dataset,distribmodelopt,tmpinverseopt);
%PrepareDataModelResidual_lola_AndPlot   (dataset,modelopt,inverseopt);
%%%%%%%%%%%%%%%%%%%%%%%%
