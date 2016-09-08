function [Lap,Lap_inv]=MakeSmoothMatrix(slip, patches, nhorz, nvert, ExtendSize)
%MakeFaultPatches             - Generates smoothing matrix (Lapacian) for distributed slip inversion
%
%  usage:  [modelopt]=MakeFaultPatches(distribopt,modelopt)
%
%  FA 7/2008 based on Sjonni Jonsson's code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %logmessage(sprintf('[]=%s(%s)',mfilename,inputname(1)));

  if size(patches,1) ==1
     Lap=0 ; Lap_inv=0;                 % case that only one patch
  else
     [Lap,Lap_inv]=modelwt2(nvert, nhorz, patches(:,1),ExtendSize(2)/nvert,1);
  end
