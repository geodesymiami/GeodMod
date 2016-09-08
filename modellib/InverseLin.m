function [modelopt]=InverseLin(dataset,inverseopt,modelopt);
%function [dislocations,parplane,parinfls]=DistOp(dg1,dataset,distopopt,parinfl);
% InverseLin-  linear inversion for fault slip and opening 
%
% usage: [mhat,mstar]=InverseLin(dataset,inverseopt,objfuncopt)
% usage: [dislocations,parplane,parinfls]=DistOp(dg1,dataset,distopopt)
%
%Input:   dataset       structure containing input data              
%         inverseopt    structure with options for DistOp
% Input:   dg1           dislocation geometry
%          dataset       structure containing input data
%          distopopt     structure with options for DistOp
%          parinfl       inflation source geometry
%
% Output:  dislocations  dislocations with opening
%          parplane      Parameter of planes
%          parinfls      Inflation source geometry and strength
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
% FA, May 2005, based on code from SJ for the DSF earthquake and for the 2 Iceland quakes
% V1.0  Falk Amelung, June 2007   based on Sjonni Jonnson's model_step1,2,3 programs
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
        'DoInverseLin'      ,   'off'                ,     ...
        'FactorLin'         ,   'off'                ,     ...
        'FactorNonLin'      ,   'off'                ,     ...
        'FactorNonLinDelta' ,    0.2                 ,     ...
        'PhaseRamp'         ,    'on'                ,     ...
        'ExtendToSurface'   ,   'on'                 ,     ...
        'PatchSize'         ,  [ 1.0  1.0]           ,     ...
        'ExtendSize'        ,   'off'                ,     ...
        'kappa'             ,    0.1                 ,     ...
        'Inflation'         ,   'off'                ,     ...
        'Factor'            ,   'off'                ,     ...
        'InverseSign'       ,   'off'                ,     ...
        'Plot'              ,   'off'      )         ;

[inverseopt]=process_defaultoptions(inverseopt,defaultopt);  
f=fieldnames(inverseopt) ; for i=1:length(f) eval([char(f{i}) '= inverseopt.(f{i}) ;' ]) ; end

if  ~DoInverseLin  return; end

if ~isfield(modelopt.par,'xy') || length(modelopt.par.xy)< 10
   logmessage('no dislocation, no inversion for slip/opening distribution possible')
   return
end

out_name            = fullfile(dir_out,'InverseLin'); 
inverseopt.out_name = out_name ;
nu                  = modelopt.dislocopt.nu;      % Poisson's ratio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dataset]      = InitializeDataset(dataset,inverseopt) ;                                       %adds G_phaseramp and normalization to dataset
[d,coord,normalization,radarlook,datind,hgt,G_phaseramp,D_1,D_2,D_3,D_4,D_5,D_6,D_7,D_8]  = datasetstructure2data(dataset);
%readfrom_dataset_structure;      % need to get rid of the following commands
 sigphi_1 = dataset(1).sigphi;sigphi_2= dataset(2).sigphi; sigphi_3 = dataset(3).sigphi;
 Ndata_1  = dataset(1).Ndata ; Ndata_2 = dataset(2).Ndata ; Ndata_3 = dataset(3).Ndata ;
 DataSet_1=dataset(1).DataSet;DataSet_2=dataset(2).DataSet;DataSet_3=dataset(3).DataSet; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup fault geometry (build patches from iput dislocations)
%

   dg1=modelopt.par.xy(1:10);

   if ExtendToSurface     dg1 = Extend2Surface(dg1);  end    % Extend Input dislocation to surface

   if ~PatchSize          PatchSize =[dg1(1) dg1(2)]; end    %set to length,width of original dislocation if 'off'
   if ~ExtendSize         ExtendSize=[dg1(1) dg1(2)]; end    %set to length,width of original dislocation if 'off'

   [lenpatch,widpatch] = deal(PatchSize(1), PatchSize(2));
   dg1(1:2)            = ExtendSize;

   len1     = dg1(1);          
   wid1     = dg1(2);
   nhE1     = len1/lenpatch;   
   nvE1     = wid1/widpatch;
   nhe1     = round(nhE1);     
   nve1     = round(nvE1);
   nel1     = nhe1.*nve1;      
   Nfaultp1 = length(nhe1);
   nel      = nel1 ;

   pm1 = patchfault(dg1(1:7)',nhe1,nve1);
   pm  = [pm1];

%% Construct kernel and data vector (weighted)
%  unweighted combined kernel

% We'll calculate opening kernel

  slip = [0,0,1];
  
  [K1,K2,K3] = MakeRadarKern(pm1, slip, coord,   nu, nhe1, nve1, radarlook', length(coord)); %radarlook needs dimension of coord
 
% Construct the kernel (or design matrix):

  Kuw = [K3];                       
  dsq = [0 0 1];                               % No dip-slip
  K   = Kuw .* repmat(1./sqrt(normalization),1,size(Kuw,2))

% Define the datavector and weight the data

 % d = [datavec_1'; datavec_2' ; datavec_3'; datavec_4' ; datavec_5'];
  if InverseSign  d=-1*d;  end                      % changing sign of data 
  dprime = d./sqrt(normalization);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct smoothing matrix
%  Create second order finite difference matrix for smoothing (laplacian)
%
  if size(pm1,1) ==1
     Lap1=0 ; Lap_inv1=0;                 % case that only one patch
  else
     [Lap1,Lap_inv1]=modelwt2(nve1, nhe1, pm1(:,1),dg1(2)/nve1,1);
  end
  
  Lap  = Lap1;
  kLap = kappa*Lap1;

% check for dip-slip

  if dsq(1) == 1;  Lap1t = BlockDiag(Lap,Lap1); Lap = Lap1t; Lap1t = BlockDiag(kLap,kappa1*Lap1); kLap = Lap1t; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% testing new geodmod code
[Gout]=MakeDesignMatrix(dataset,[],inverseopt);  %G=Gout;
Gtop = [ K    Gout];
Gsmo = [ kLap zeros(size(Lap,1),size(Gout,2)) ];
G    = [Gtop  ; Gsmo];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solve using fast non-negative least squares

din  = [dprime; zeros(size(K,2),1)];


sal = fnnls(G'*G,G'*din);


s     = sal;  % the slip of all faults

%% RMS calculations (data prediction and residual)
%
pred =  K*s ./sqrt(normalization);

resi = dprime-pred;
resi = resi.*sqrt(normalization);

%sqrflag=true;
%[resi,pred,u,rms,npar,weights,rms_unitsig,mlin,dataset] = feval(objfunc,par,dataset,objfuncopt,sqrflag);                     % FA Feb 2008. Needed to make it work for ProjectProfile


    tot_rms  = norm(resi)                            / sqrt(length(resi));
    rms_1    = norm(resi(1:datind(1)))                 / sqrt(Ndata_1);
if D_2 rms_2 = norm(resi(datind(1)+1:datind(2))) / sqrt(Ndata_2); end
if D_3 rms_3 = norm(resi(datind(2)+1:datind(3))) / sqrt(Ndata_3); end
if D_4 rms_4 = norm(resi(datind(3)+1:datind(4))) / sqrt(Ndata_4); end
if D_5 rms_5 = norm(resi(datind(4)+1:datind(5))) / sqrt(Ndata_5); end

        St = sprintf('RMS              : %2.2f cm',tot_rms*100); disp(St)
        St = sprintf('RMS dataset %s   : %2.2f cm',DataSet_1,rms_1*100); disp(St)
if D_2  St = sprintf('RMS dataset %s   : %2.2f cm',DataSet_2,rms_2*100); disp(St) ; end
if D_3  St = sprintf('RMS dataset %s   : %2.2f cm',DataSet_3,rms_3*100); disp(St) ; end
if D_4  St = sprintf('RMS dataset %s   : %2.2f cm',DataSet_4,rms_4*100); disp(St) ; end
if D_5  St = sprintf('RMS dataset %s   : %2.2f cm',DataSet_5,rms_5*100); disp(St) ; end

% Calculate the weight of each data set
%Wper = [sum(W(1:Ndata_1)) sum(W(Ndata_1+1:Ndata_1+Ndata_2))]';
        W    = 1./sqrt(normalization);
        Wper = [       sum(W(          1:datind(1)))];
if D_2  Wper = [Wper ; sum(W(datind(1)+1:datind(2)))];  end ;
if D_3  Wper = [Wper ; sum(W(datind(2)+1:datind(3)))];  end ;
if D_4  Wper = [Wper ; sum(W(datind(3)+1:datind(4)))];  end ;
if D_5  Wper = [Wper ; sum(W(datind(4)+1:datind(5)))];  end ;
        Wper = Wper/sum(Wper)*100;
disp(' ')
       St = sprintf('Weight dataset %s: %2.1f%s',DataSet_1,Wper(1),'%'); disp(St)
if D_2 St = sprintf('Weight dataset %s: %2.1f%s',DataSet_2,Wper(2),'%'); disp(St); end
if D_3 St = sprintf('Weight dataset %s: %2.1f%s',DataSet_3,Wper(3),'%'); disp(St); end
if D_4 St = sprintf('Weight dataset %s: %2.1f%s',DataSet_4,Wper(4),'%'); disp(St); end
if D_5 St = sprintf('Weight dataset %s: %2.1f%s',DataSet_5,Wper(5),'%'); disp(St); end

% Calculate the solution roughness
tmp1=Lap*s;
rough=sum(abs(tmp1))/(2*(nel));

disp(' ')
St = sprintf('%s%1.2f%s','Mean roughness   : ',rough*100,' cm/km'); disp(St);
%% compute the moment of the estimated model

% number of fault patches
nf=size(pm,1);

% pick out strike slip from solution
sslip = s(1:nf);

if slip(1)==1
   dslip=s(nel+1:2*nel);
else
   dslip=zeros(size(sslip));
end
	
% Define the rigidity
rig = 3.0e10;

% Calculate slip magnitude	
   Nels = nel;
   if length(s) >= 2*Nels
      sslipmag = abs(s(1:Nels));
      dslipmag = abs(s(Nels+1:2*Nels));
      slipmag  = ( s(1:Nels).^2 + s(Nels+1:2*Nels).^2 ).^(0.5);
   else
      slipmag = abs(s(1:Nels));
   end
   
   
% calculate "geometric" moment of dip slip and strike slip components    
   for j = 1:Nels
     if length(s) >= 2*Nels
       ssgmom(j) = pm(j,1)*pm(j,2)*1.0e6*sslipmag(j);
       dsgmom(j) = pm(j,1)*pm(j,2)*1.0e6*dslipmag(j);       
       gmom(j)   = pm(j,1)*pm(j,2)*1.0e6*slipmag(j); 
     else
       gmom(j) = pm(j,1)*pm(j,2)*1.0e6*slipmag(j);
     end	%convert to meters squared
   end
   
% Calculate total moment and moment magnitude
   %ind1 = [1:nel1]; ind2 = [nel1+1:nel1+nel2];  
   moment = rig.*[sum(gmom) ];
   Mw = 2*log10(moment)/3 - 6.0;
   totalop=sum(gmom)*1.0e-6 ;
 
% Calculate the ratio between dip-slip and strike-slip moment,
% and print on the sceen
   if length(s) >= 2*Nels
     ssmom1 = rig*sum(ssgmom); 
     dsmom1 = rig*sum(dsgmom); 
     ssper1 = ssmom1/(ssmom1+dsmom1)*100; 
     dsper1 = dsmom1/(ssmom1+dsmom1)*100; 
     msl1 = max(dslip); 
     S = sprintf('%s%1.1f%s','Strike-slip moment       : ',ssper1,' %'); disp(S); 
     S = sprintf('%s%1.1f%s','Dip-slip moment          : ',dsper1,' %'); disp(S);
     disp(' ')        
     S = sprintf('%s%1.2f%s','Maximum dip-slip         : ',msl1,' m'); disp(S);
   end
   
   msl1 = max(sslip); 

   if slip(3)==1;
      S = sprintf('%s%1.2f%s','Maximum opening  : ',msl1,' m'); disp(S);
      S = sprintf('%s%1.3g%s','Total opening    : ',totalop,' 10^6 m^3'); disp(S);
   else
      S = sprintf('%s%1.2f%s','Maximum Strike-slip      : ',msl1,' m'); disp(S);
      S = sprintf('%s%1.3g%s',    'Total seismic moment  Mo = ',moment(1),' Nm'); disp(S);
      S = sprintf('%s%1.2f','Moment magnitude      Mw = ',Mw(1)); disp(S);
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   parinfls=[];
   if InverseSign s=-s; parplane=-parplane;  end
   if InverseSign && inflation infl_strength=-infl_strength ;  end
     %if inflation   parinfls=[pinfl' ; infl_strength] ; end

   tmp_dislocations=zeros(size(pm,1),10) ;
   tmp_dislocations(1:length(pm(:)))=pm(:);
   tmp_dislocations(:,10)=s;
   tmp_dislocations=tmp_dislocations';
   dislocations=tmp_dislocations(:);

plot=true;
if ~plot return ; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------
% Script to plot up the estimated slip distribution as patches in 3D
%-------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load colormap
load slipcol

%plot fault geometry
[fx,fy,fz]=flakes(pm');
figure; fill3(fx,fy,-fz,[1:size(pm,1)]); axis image

%--------------------------------------------------------
% if Dip-slip plot separetely
%if length(s)>nf+12
if slip(1)==1      
      % Plot Strike slip and dip slip distribution
      %dslip = s(nf+1:2*nf);
      figure
      subplot(211)
      fill3(fx,fy,-fz,sslip'); colormap(slipcol); colorbar; axis image
      title(['Predicted Strike-slip Distribution, roughness = ', num2str(100*rough),' cm/km']); hold on
      %plot3(hypo(:,1),hypo(:,2),hypo(:,3),'*k')	

      subplot(212)
      fill3(fx,fy,-fz,dslip'); colormap(slipcol); colorbar; axis image
      title(['Predicted Dip-slip Distribution, roughness = ', num2str(100*rough),' cm/km']); hold on
      %plot3(hypo(:,1),hypo(:,2),hypo(:,3),'*k')	
      
      % Plot Slip magnitude and Rake

      slipmag  = (sslip.^2 + dslip.^2).^(1/2);
      sliprake = 180 - atan(dslip./(sslip+1e-6))*180/pi;
      
      figure
      subplot(211)
      fill3(fx,fy,-fz,slipmag'); colormap(slipcol); colorbar; axis image
      title('Predicted Slip Magnitude'); hold on
      %plot3(hypo(:,1),hypo(:,2),hypo(:,3),'*k')	

      subplot(212)
      fill3(fx,fy,-fz,sliprake'); colormap(slipcol); colorbar; axis image
      title('Predicted Rake'); hold on
      %plot3(hypo(:,1),hypo(:,2),hypo(:,3),'*k')	

    
else  % IF strike slip only
      figure      
      fill3(fx,fy,-fz,sslip'); colormap(slipcol); colorbar; axis image;
      title(['Predicted Strike-slip Distribution, Roughness = ',num2str(rough*100), ' cm/km']); hold on
      %plot3(hypo(:,1),hypo(:,2),hypo(:,3),'*k')	
end


nn=[nel1]; ca = [0 max(abs(s))];
if slip(3)==1;
   [ffx,ffz]=flakes2D(pm',nhe1,nve1);
   figure;set(gcf,'Position',[428 68 460 560]); colormap(slipcol)
   subplot(211); patch(ffx,ffz,sslip'); axis image; hold on; caxis(ca);
   title('North     (strike-slip)   South'); colorbar;ax=axis;%plot(hypo(1,2),hypo(1,3),'*');ax=axis;
   ylabel('Depth down-dip (km)'); xlabel('Distance along fault')

   title('South     (opening)   North'); colorbar;ax=axis;%plot(hypo(1,2),hypo(1,3),'*');ax=axis;
   ylabel('Depth (km)'); xlabel('Distance along riftzone')
   subplot(212);  patch(ffx,ffz,dslip'); hold on
   axis image; title('North   (dip-slip)    South'); caxis(ca); colorbar;%plot(hypo(2,2),hypo(2,3),'*')   
   ylabel('Depth down-dip (km)'); xlabel('Distance along fault')
end
