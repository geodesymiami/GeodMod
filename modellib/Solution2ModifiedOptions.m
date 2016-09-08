function [modinverseopt,modmodelopt] = Solution2ModifiedOptions(N_disloc,modelpar,inverseopt,dataset,roughness)
%   GenerateSummary       - Calculates RMS, Weights, and writes them into a string                  
%   FA 7/2008
%

  
  modmodelopt.Unit               = inverseopt.objfuncopt.modelopt.Unit;
  modmodelopt.N_disloc           = N_disloc;


  modinverseopt.algorithm        = 'Distributed Slip';
  modinverseopt.objfunc          = inverseopt.objfunc;
  modinverseopt.FollowGradient   = inverseopt.FollowGradient;
  modinverseopt.plotdataopt      = inverseopt.plotdataopt;
  modinverseopt.objfuncopt       = inverseopt.objfuncopt;
  modinverseopt.distribopt       = inverseopt.distribopt;
  modinverseopt.FactorLin        = inverseopt.FactorLin;
  modinverseopt.PhaseRamp        = inverseopt.PhaseRamp;
  modinverseopt.FactorNonLin     = inverseopt.FactorNonLin;
  modinverseopt.FactorNonLinDelta= inverseopt.FactorNonLinDelta;
  modinverseopt.plotmodelopt     = inverseopt.plotmodelopt;
  modinverseopt.disloc_bounds.xy = [modelpar(:)  modelpar(:)];
  modmodelopt.par.xy=modelpar(:);
  [modmodelopt]                  = InitializeModelopt(modmodelopt,inverseopt.plotdataopt.basemap);
  [qmodinverseopt]               = InitializeInverseopt(dataset,modinverseopt,modmodelopt); 

  modinverseopt = qmodinverseopt;

  if exist('roughness','var') modmodelopt.roughness = roughness; end
