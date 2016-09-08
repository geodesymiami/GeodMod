function [phaseramp] = phaseramp_dummy(dataset,objfuncopt)
%   phaseramp_dummy  - generates a dummy phaseramp to remove from the data so that fnnls works in InverseLinDistrib   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [SAR_1,SAR_2,SAR_3,SAR_4,SAR_5,SAR_6,SAR_7,SAR_8]  = dealfalk(dataset.SAR);
  PhaseRamp                                          = objfuncopt.PhaseRamp ;
  
% calculate x and y coordinates for plane removal from data_grid

switch PhaseRamp
case{'Ramp',true}
   dummy = [-1 -0.01 -0.01]' ;
case{'Const'}
   dummy = -1;
otherwise
   dummy =[];
end

phaseramp = repmat(dummy,sum([dataset.SAR]), 1);
