function [dataset] = ModifyDatasetLin(dataset,objfuncopt,m_phaseramp,replace)
%   ModifyDatasetLin  - Modify Dataset using  linear model parameters (plane,const parameters)  
%
% INPUT:
%       mlin                  - linear model parameters (first source strength, then plane,const parameters 
%       objfunctopt           - options structure with following fields:
%            PhaseRamp          'plane' or 'const'
%            FactorNonLin(??)   'on' if invert for a factor for each interferogram
%
%  used in GenericObjectiveFunction (only)
%
%  FA June 2007. Verified that it works fine for PhaseRamp option
%                The FactorNonLin or SARmul,GPSmul is not incorporated here.
%                It never has verified
%  FA June 2007. Now puts modified data into dataset.data_mod
%  FA Feb 2008 comment: if mlin empty we simply should have dataset_data_mod=dataset.data; 
%                       (currently we have mlin=1 even if mogi strength is fixed)
%  FA September 2009: added 'replace' option
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if exist('replace','var') && strcmp(replace,'replace') replace_flag=true ; else replace_flag=false; end
  
  [SAR_1,SAR_2,SAR_3,SAR_4,SAR_5,SAR_6,SAR_7,SAR_8]  = dealfalk(dataset.SAR);
  
  if ~SAR_1 return; end  %FA 11/2008 if only GPS 
  
  PhaseRamp = objfuncopt.PhaseRamp ;
  FactorLin = objfuncopt.FactorLin ;
  
% calculate x and y coordiantes for plane removal from data_grid

if SAR_1
   data_1=dataset(1).data; x_posting   =dataset(1).x_posting;  y_posting   =dataset(1).y_posting;
   x_l=([1:size(data_1,2)] - 1)*x_posting;   %axis label vector
   y_l=([1:size(data_1,1)] - 1)*y_posting;   %axis label vector
   [xx,yy]=meshgrid(x_l,y_l); planegrid = [ones(length(xx(:)),1)  xx(:) yy(:)] ;
   [dataset(:).data_mod]=deal(dataset(:).data);        % initialize so that *.data_mod can be filled
end

switch PhaseRamp
case{'Ramp',true}
   if SAR_1 dataset(1).data_mod(:) = dataset(1).data(:)-(planegrid*m_phaseramp( 1: 3));  end
   if SAR_2 dataset(2).data_mod(:) = dataset(2).data(:)-(planegrid*m_phaseramp( 4: 6));  end 
   if SAR_3 dataset(3).data_mod(:) = dataset(3).data(:)-(planegrid*m_phaseramp( 7: 9));  end 
   if SAR_4 dataset(4).data_mod(:) = dataset(4).data(:)-(planegrid*m_phaseramp(10:12));  end 
   if SAR_5 dataset(5).data_mod(:) = dataset(5).data(:)-(planegrid*m_phaseramp(13:15));  end 
  
   if SAR_1 dataset(1).datavec_mod = dataset(1).datavec-([ones(size(dataset(1).coord,2),1) dataset(1).coord']*m_phaseramp( 1: 3))'; end
   if SAR_2 dataset(2).datavec_mod = dataset(2).datavec-([ones(size(dataset(2).coord,2),1) dataset(2).coord']*m_phaseramp( 4: 6))'; end
   if SAR_3 dataset(3).datavec_mod = dataset(3).datavec-([ones(size(dataset(3).coord,2),1) dataset(3).coord']*m_phaseramp( 7: 9))'; end
   if SAR_4 dataset(4).datavec_mod = dataset(4).datavec-([ones(size(dataset(4).coord,2),1) dataset(4).coord']*m_phaseramp(10:12))'; end
   if SAR_5 dataset(5).datavec_mod = dataset(5).datavec-([ones(size(dataset(5).coord,2),1) dataset(5).coord']*m_phaseramp(13:15))'; end
case{'Const'}
   if SAR_1 dataset(1).data_mod(:) = dataset(1).data(:)-m_phaseramp(1);  end
   if SAR_2 dataset(2).data_mod(:) = dataset(2).data(:)-m_phaseramp(2);  end
   if SAR_3 dataset(3).data_mod(:) = dataset(3).data(:)-m_phaseramp(3);  end
   if SAR_4 dataset(4).data_mod(:) = dataset(4).data(:)-m_phaseramp(4);  end
   if SAR_5 dataset(5).data_mod(:) = dataset(5).data(:)-m_phaseramp(5);  end
    
   if SAR_1 dataset(1).datavec_mod = dataset(1).datavec-m_phaseramp(1); end
   if SAR_2 dataset(2).datavec_mod = dataset(2).datavec-m_phaseramp(2); end
   if SAR_3 dataset(3).datavec_mod = dataset(3).datavec-m_phaseramp(3); end
   if SAR_4 dataset(4).datavec_mod = dataset(4).datavec-m_phaseramp(4); end
   if SAR_5 dataset(5).datavec_mod = dataset(5).datavec-m_phaseramp(5); end
otherwise
   if SAR_1 dataset(1).data_mod(:) = dataset(1).data(:); end
   if SAR_2 dataset(2).data_mod(:) = dataset(2).data(:); end
   if SAR_3 dataset(3).data_mod(:) = dataset(3).data(:); end
   if SAR_4 dataset(4).data_mod(:) = dataset(4).data(:); end
   if SAR_5 dataset(5).data_mod(:) = dataset(5).data(:); end
       
   if SAR_1 dataset(1).datavec_mod = dataset(1).datavec; end
   if SAR_2 dataset(2).datavec_mod = dataset(2).datavec; end
   if SAR_3 dataset(3).datavec_mod = dataset(3).datavec; end
   if SAR_4 dataset(4).datavec_mod = dataset(4).datavec; end
   if SAR_5 dataset(5).datavec_mod = dataset(5).datavec; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                                                
%% NOTE: For fac not possible because it applys to the  %%                                                                                                                
%% sources (Disloc and mogi ) and not to the datasets   %%                                                                                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                                                
                                                                                                                                                                          
%[GPS_1,GPS_2,GPS_3,GPS_4,GPS_5,GPS_6,GPS_7,GPS_8]  = dealfalk(dataset.GPS);
%[ind_1,ind_2,ind_3,ind_4,ind_5,ind_6,ind_7,ind_8] = dealfalk(dataset.ind);        

if FactorLin                                                                                                                                                              
   %if GPS    Gf = blkdiag(Gf, G(ind_GPS,:))*0           ;        
   %           if SAR_1  Gf = blkdiag(Gf, G(ind_1  ,:))   ;  end
   % else                Gf = blkdiag(Gf, G(ind_1  ,:))*0 ;  end
   %                                                                                                                                                                      
   % if SAR_2  Gf = blkdiag(Gf        , G(ind_2,:))       ; end
   % if SAR_3  Gf = blkdiag(Gf        , G(ind_3,:))       ; end
   % if SAR_4  Gf = blkdiag(Gf        , G(ind_4,:))       ; end
   % if SAR_5  Gf = blkdiag(Gf        , G(ind_5,:))       ; end
end  

% June 2007: The following does not take data_mod into account. Should be easy to incorporate (in datasetstructure2data.m check use                                       
% d = dataset(:).datavec_mod if 'DataField','data_mod'). Have not done this because fulldata it does not look used to me for plotting 

if replace_flag
   if SAR_1  i=1; dataset(i).data(:)=dataset(i).data_mod(:);  dataset(i).datavec=dataset(i).datavec_mod; end
   if SAR_2  i=2; dataset(i).data(:)=dataset(i).data_mod(:);  dataset(i).datavec=dataset(i).datavec_mod; end
   if SAR_3  i=3; dataset(i).data(:)=dataset(i).data_mod(:);  dataset(i).datavec=dataset(i).datavec_mod; end
   if SAR_4  i=4; dataset(i).data(:)=dataset(i).data_mod(:);  dataset(i).datavec=dataset(i).datavec_mod; end
   if SAR_5  i=5; dataset(i).data(:)=dataset(i).data_mod(:);  dataset(i).datavec=dataset(i).datavec_mod; end
   dataset = rmfield(dataset,{'data_mod' 'datavec_mod'});
   if isfield(dataset,'fulldata') dataset= rmfield(dataset,'fulldata'); end
   
end

if isfield(dataset,'fulldata')
   dataset=rmfield(dataset,'fulldata');
   [d.d,d.coord,d.normalization,d.radarlook,d.datind,d.hgt,d.G_phaseramp,d.D_1,d.D_2,d.D_3,d.D_4,d.D_5,d.D_6,d.D_7,d.D_8]  =datasetstructure2data(dataset,objfuncopt);
   dataset(1).fulldata=d;
end    
