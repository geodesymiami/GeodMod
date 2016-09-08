function [Gout] = MakeDesignMatrix(dataset,G,inverseopt)
%   MakeDesignMatrix  - Design matrix for linear inversion  of slip (opening), plane and factor for each data set 
%
% INPUT:
%		G	             - existing design matrix  
%       inversetopt      - options structure with following fields:
%       
%       linear_inversion -  'plane' or 'const'
%       FactorLin   -   'on' to invert for a factor for each source (Disloc and Inflation) in each interferogram
%
% OUTPUT:
%
%

FactorLin=inverseopt.FactorLin ;
PhaseRamp=inverseopt.PhaseRamp ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[D_1,D_2,D_3,D_4,D_5,D_6,D_7,D_8]                                 = dealfalk(dataset.exist);
[coord_1,coord_2,coord_3,coord_4,coord_5,coord_6,coord_7,coord_8] = dealfalk(dataset.coord);
[Ndata_1,Ndata_2,Ndata_3,Ndata_4,Ndata_5,Ndata_6,Ndata_7,Ndata_8] = dealfalk(dataset.Ndata);
[SAR_1,SAR_2,SAR_3,SAR_4,SAR_5,SAR_6,SAR_7,SAR_8]                 = dealfalk(dataset.SAR);
[GPS_1,GPS_2,GPS_3,GPS_4,GPS_5,GPS_6,GPS_7,GPS_8]                 = dealfalk(dataset.GPS);
[ind_1,ind_2,ind_3,ind_4,ind_5,ind_6,ind_7,ind_8]                 = dealfalk(dataset.ind);

DataSet = cell(length(dataset),1);[DataSet{:}]=deal(dataset.DataSet);ieast=strmatch('GPSeast',DataSet);iup=strmatch('GPSup  ' ,DataSet);
ind_GPS = [dataset(ieast).ind dataset(ieast+1).ind dataset(iup).ind]; GPS=~isempty(ind_GPS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
[Gp,Gf]=deal([]);

switch PhaseRamp
case {'Ramp',true}
   if SAR_1  plane_1 = [ ones(Ndata_1,1)  coord_1(1:2,:)'];   Gp = blkdiag(Gp,plane_1) ; end
   if SAR_2  plane_2 = [ ones(Ndata_2,1)  coord_2(1:2,:)'];   Gp = blkdiag(Gp,plane_2) ; end
   if SAR_3  plane_3 = [ ones(Ndata_3,1)  coord_3(1:2,:)'];   Gp = blkdiag(Gp,plane_3) ; end
   if SAR_4  plane_4 = [ ones(Ndata_4,1)  coord_4(1:2,:)'];   Gp = blkdiag(Gp,plane_4) ; end
   if SAR_5  plane_5 = [ ones(Ndata_5,1)  coord_5(1:2,:)'];   Gp = blkdiag(Gp,plane_5) ; end
   if GPS    Gp      = [ Gp;zeros(length(ind_GPS),size(Gp,2))]                         ; end  %fill matrix with zeros because there is nothing to invert for
case {'Const'}
   if SAR_1  Gp = blkdiag( Gp , ones(Ndata_1,1)) ;  end 
   if SAR_2  Gp = blkdiag( Gp , ones(Ndata_2,1)) ;  end 
   if SAR_3  Gp = blkdiag( Gp , ones(Ndata_3,1)) ;  end 
   if SAR_4  Gp = blkdiag( Gp , ones(Ndata_4,1)) ;  end 
   if SAR_5  Gp = blkdiag( Gp , ones(Ndata_5,1)) ;  end 
   
   if GPS_1  Gp = blkdiag( Gp , ones(Ndata_1,1)) ;  end 
   if GPS_2  Gp = blkdiag( Gp , ones(Ndata_2,1)) ;  end 
   if GPS_3  Gp = blkdiag( Gp , ones(Ndata_3,1)) ;  end 
   if GPS_4  Gp = blkdiag( Gp , ones(Ndata_4,1)) ;  end 
   if GPS_5  Gp = blkdiag( Gp , ones(Ndata_5,1)) ;  end 
   if GPS_6  Gp = blkdiag( Gp , ones(Ndata_6,1)) ;  end 
   if GPS_7  Gp = blkdiag( Gp , ones(Ndata_7,1)) ;  end 
   if GPS_8  Gp = blkdiag( Gp , ones(Ndata_8,1)) ;  end 
   
   %if GPS    Gp = [Gp; zeros(length(ind_GPS),size(Gp,2))]; end      % 8/08:  empty if only GPS. Is that correct ?
   %if GPS    Gp = [Gp; ones(length(ind_GPS),1)]; end            % 8/08: test to invert for constant GPS offset
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

Gout       = [G Gp Gf ]  ;
