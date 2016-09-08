function [str] = MakeStringLinearInversionInfo(inverseopt,dataset,mlin)
%   MakeStringLinearInversionInfo  - Generates a string with information on linear inversion parameters
%
% INPUT:
%       objfuncopt       - options structure for objective function
%       mlin             - parameter with linear inversion parameters (obtained by calling objective function)  
%
% OUTPUT:
%       str              - string with info on linear inversion parameters
%
% FA: Feb 2008
% TODO: This function is very similar to MakeDesignMatrix and the functionality should be merged (performance?)

FactorLin   = inverseopt.FactorLin ;
PhaseRamp   = inverseopt.PhaseRamp ;
linearind   = inverseopt.linearind ;
objfuncopt  = inverseopt.objfuncopt;
m_phaseramp = mlin(length(linearind)+1:end);         % linear phaseramp  parameters

[SAR_1,SAR_2,SAR_3,SAR_4,SAR_5,SAR_6,SAR_7,SAR_8]=dealfalk(dataset.SAR);
[GPS_1,GPS_2,GPS_3,GPS_4,GPS_5,GPS_6,GPS_7,GPS_8]=dealfalk(dataset.GPS);

DataSet=cell(length(dataset),1);[DataSet{:}]=deal(dataset.DataSet);ieast=strmatch('GPSeast',DataSet);iup=strmatch('GPSup  ',DataSet);
ind_GPS=[dataset(ieast).ind dataset(ieast+1).ind dataset(iup).ind]; GPS=~isempty(ind_GPS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


str = sprintf('Number of model parameters: %d linear, %d fixed, %d free, %d phaseramp\n',length(objfuncopt.linearind),length(objfuncopt.fixind),length(objfuncopt.freeind), length(m_phaseramp));

switch PhaseRamp
case {'Ramp',true}
   if SAR_1  str = [str sprintf('%-7s Ramp: %6.4f %6.4f %6.4f\n' ,dataset(1).DataSet,m_phaseramp(1:3))]; m_phaseramp(1:3)=[]; end
   if SAR_2  str = [str sprintf('%-7s Ramp: %6.4f %6.4f %6.4f\n' ,dataset(2).DataSet,m_phaseramp(1:3))]; m_phaseramp(1:3)=[]; end
   if SAR_3  str = [str sprintf('%-7s Ramp: %6.4f %6.4f %6.4f\n' ,dataset(3).DataSet,m_phaseramp(1:3))]; m_phaseramp(1:3)=[]; end
   if SAR_4  str = [str sprintf('%-7s Ramp: %6.4f %6.4f %6.4f\n' ,dataset(4).DataSet,m_phaseramp(1:3))]; m_phaseramp(1:3)=[]; end
   if SAR_5  str = [str sprintf('%-7s Ramp: %6.4f %6.4f %6.4f\n' ,dataset(5).DataSet,m_phaseramp(1:3))]; m_phaseramp(1:3)=[]; end
   if GPS  ; end             % no plane
case {'Const'}
   if SAR_1  str = [str sprintf('%-7s Const: %6.4f\n' ,dataset(1).DataSet,m_phaseramp(1))]; m_phaseramp(1)=[]; end
   if SAR_2  str = [str sprintf('%-7s Const: %6.4f\n' ,dataset(2).DataSet,m_phaseramp(1))]; m_phaseramp(1)=[]; end
   if SAR_3  str = [str sprintf('%-7s Const: %6.4f\n' ,dataset(3).DataSet,m_phaseramp(1))]; m_phaseramp(1)=[]; end
   if SAR_4  str = [str sprintf('%-7s Const: %6.4f\n' ,dataset(4).DataSet,m_phaseramp(1))]; m_phaseramp(1)=[]; end
   if SAR_5  str = [str sprintf('%-7s Const: %6.4f\n' ,dataset(5).DataSet,m_phaseramp(1))]; m_phaseramp(1)=[]; end
   
   if GPS_1  str = [str sprintf('%-7s Const: %6.4f\n' ,dataset(1).DataSet,m_phaseramp(1))]; m_phaseramp(1)=[]; end  % 8/08: Testing for GPS const
   if GPS_2  str = [str sprintf('%-7s Const: %6.4f\n' ,dataset(2).DataSet,m_phaseramp(1))]; m_phaseramp(1)=[]; end
   if GPS_3  str = [str sprintf('%-7s Const: %6.4f\n' ,dataset(3).DataSet,m_phaseramp(1))]; m_phaseramp(1)=[]; end
   if GPS_4  str = [str sprintf('%-7s Const: %6.4f\n' ,dataset(4).DataSet,m_phaseramp(1))]; m_phaseramp(1)=[]; end
   if GPS_5  str = [str sprintf('%-7s Const: %6.4f\n' ,dataset(5).DataSet,m_phaseramp(1))]; m_phaseramp(1)=[]; end
   if GPS_6  str = [str sprintf('%-7s Const: %6.4f\n' ,dataset(6).DataSet,m_phaseramp(1))]; m_phaseramp(1)=[]; end
   if GPS_7  str = [str sprintf('%-7s Const: %6.4f\n' ,dataset(7).DataSet,m_phaseramp(1))]; m_phaseramp(1)=[]; end
   if GPS_8  str = [str sprintf('%-7s Const: %6.4f\n' ,dataset(8).DataSet,m_phaseramp(1))]; m_phaseramp(1)=[]; end
   if GPS  ; end             % no constant
end

if FactorLin
   if GPS    Gf = blkdiag(Gf, G(ind_GPS,:))*0           ; 
             if SAR_1  Gf = blkdiag(Gf, G(ind_1  ,:))   ;  end
   else                Gf = blkdiag(Gf, G(ind_1  ,:))*0 ;  end

   if SAR_2  Gf = blkdiag(Gf        , G(ind_2,:)); end
   if SAR_3  Gf = blkdiag(Gf        , G(ind_3,:)); end
   if SAR_4  Gf = blkdiag(Gf        , G(ind_4,:)); end
   if SAR_5  Gf = blkdiag(Gf        , G(ind_5,:)); end
end

