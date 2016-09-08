function [dataset] = ModifyDatasetNonLin(par,dataset,objfuncopt,mlin)
%   ModifyDataseNontLin  - Modify Dataset according to results of linear inversion  
%
% INPUT:
%		G	             - existing design matrix  
%       objfunctopt      - options structure with following fields:
%       
%       linear_inversion -  'plane' or 'const'
%       factor           -   'on' if invert for a factor for each interferogram
%
% OUTPUT:
%
%

FactorLin      = objfuncopt.FactorLin ;
PhaseRamp      = objfuncopt.PhaseRamp ;
StrengthFactor = strcmp('on',objfuncopt.StrengthFactor) ;

N_disloc    = floor(length(par)/10);
if N_disloc == 0 N_disloc=1 ; end     % this accounts for mlin(1)=0 if there is no dislocation (this is because pred_disloc=0*pred in ManyDislocOneInflation)

if strcmp('off',objfuncopt.inflation)  plstart=N_disloc+1;   %get index of start of first plane in mlin
   else                                plstart=N_disloc+2;
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[D_1,D_2,D_3,D_4,D_5]                    =dealfalk(dataset.exist);
[SAR_1,SAR_2,SAR_3,SAR_4,SAR_5]          =dealfalk(dataset.SAR);
[coord_1,coord_2,coord_3,coord_4,coord_5]=dealfalk(dataset.coord);
[Ndata_1,Ndata_2,Ndata_3,Ndata_4,Ndata_5]=dealfalk(dataset.Ndata);
[ind_1,ind_2,ind_3,ind_4,ind_5,ind_6,ind_7,ind_8]=dealfalk(dataset.ind);
if ~SAR_1 return ; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% calculate x and y coordiantes for plane removal from data_grid
data_1=dataset(1).data; x_posting   =dataset(1).x_posting;  y_posting   =dataset(1).y_posting;
x_l=([1:size(data_1,2)] - 1)*x_posting;   %axis label vector
y_l=([1:size(data_1,1)] - 1)*y_posting;   %axis label vector
[xx,yy]=meshgrid(x_l,y_l); planegrid = [ones(length(xx(:)),1)  xx(:) yy(:)] ;

switch PhaseRamp
case{'Ramp',true}
   if SAR_1  plane = [ ones(Ndata_1,1)  coord_1(1:2,:)']; dataset(1).datavec=dataset(1).datavec-(plane*mlin(plstart   :plstart+2  ))'; end
   if SAR_2  plane = [ ones(Ndata_2,1)  coord_2(1:2,:)']; dataset(2).datavec=dataset(2).datavec-(plane*mlin(plstart+3 :plstart+5  ))'; end
   if SAR_3  plane = [ ones(Ndata_3,1)  coord_3(1:2,:)']; dataset(3).datavec=dataset(3).datavec-(plane*mlin(plstart+6 :plstart+8  ))'; end
   if SAR_4  plane = [ ones(Ndata_4,1)  coord_4(1:2,:)']; dataset(4).datavec=dataset(4).datavec-(plane*mlin(plstart+9 :plstart+11 ))'; end
   if SAR_5  plane = [ ones(Ndata_5,1)  coord_5(1:2,:)']; dataset(5).datavec=dataset(5).datavec-(plane*mlin(plstart+12:plstart+14))';  end

   if SAR_1  dataset(1).data(:)=dataset(1).data(:)-(planegrid*mlin(plstart   :plstart+2  ));  end
   if SAR_2  dataset(2).data(:)=dataset(2).data(:)-(planegrid*mlin(plstart+3 :plstart+5  ));  end
   if SAR_3  dataset(3).data(:)=dataset(3).data(:)-(planegrid*mlin(plstart+6 :plstart+8  ));  end
   if SAR_4  dataset(4).data(:)=dataset(4).data(:)-(planegrid*mlin(plstart+9 :plstart+11 ));  end
   if SAR_5  dataset(5).data(:)=dataset(5).data(:)-(planegrid*mlin(plstart+12:plstart+14));   end
case{'Const')
errordlg('not yet programmed for inverseopt.PhaseRamp='Const');error
%   if SAR_1  Gp =                                            ones(Ndata_1,1);     end
%   if SAR_2  Gp = [Gp zeros(size(Gp,1),1); zeros(Ndata_2,1)  ones(Ndata_2,1)];    end 
%   if SAR_3  Gp = [Gp zeros(size(Gp,1),1); zeros(Ndata_3,2)  ones(Ndata_3,1)];    end 
%   if SAR_4  Gp = [Gp zeros(size(Gp,1),1); zeros(Ndata_4,3)  ones(Ndata_4,1)];    end 
%   if SAR_5  Gp = [Gp zeros(size(Gp,1),1); zeros(Ndata_5,4)  ones(Ndata_5,1)];    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOTE: NEEDS TO BE PROGRAMMED FOR const             %%
%% For   fac   not possible because it applys to the  %%
%% sources (Disloc and mogi ) and not to the datasets %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if StrengthFactor

%   if GPS    Gf = blkdiag(Gf, G(ind_GPS,:))*0           ;
%             if SAR_1  Gf = blkdiag(Gf, G(ind_1  ,:))   ;  end
%   else                Gf = blkdiag(Gf, G(ind_1  ,:))*0 ;  end
%
%   if SAR_2  Gf = blkdiag(Gf        , G(ind_2,:))       ; end
%   if SAR_3  Gf = blkdiag(Gf        , G(ind_3,:))       ; end
%   if SAR_4  Gf = blkdiag(Gf        , G(ind_4,:))       ; end
%   if SAR_5  Gf = blkdiag(Gf        , G(ind_5,:))       ; end
 
end

