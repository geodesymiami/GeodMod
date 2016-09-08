function [u_radarlook,u,u_sum] = ForwardModel(par,coord,radarlook,modelopt,hgt)
%   ForwardModel   - calculates surface displacement for one or multiple sources 
%
%      usage:  [u_radarlook,u]       =  ForwardModel(par,coord,modelopt);
%              [u_radarlook,u,u_sum] =  ForwardModel(par,coord,modelopt);
%
% INPUT:
%		par	     - (m-by-1) vector of modelparameters  (e.g. 34-by-1 for 3 dislocations and 1 Mogi source)
%                  required order of sources: disloc,mogi,penny,mctigue,yang,visco1d,lockedandcreep
%		coord	 - (2-by-n) matrix of coordinates
%       modelopt - structure containing model options
%
% OUTPUT
%       u_radarlook - displacement in radarlook direction (1 column for each source, Ndata-by-Nsources)   
%		u           - column vector of surface displacement (3*Ndata-by-Nsources)
%		u_sum       - summed surface displacement (3*Ndata x 1)
%
%       Currently supported sources:
%
%       disloc     [ Len   Wid   Dep   Dip    Strike  xE      xN      ss     ds    op ]
%       mogi       [ xE    xN    Dep   Stren                                          ]
%       penny      [ xE    xN    Dep   Rad    Stren                                   ]
%       mctigue    [ xE    xN    Dep   Rad    Press                                   ]
%       yang       [ xE    xN    Dep   Press  majAx   AxRatio Strike  Plunge          ]
%       visco1d    [ Len   Wid   Dep   Dip    Strike  xE      xN      ss      ds   op ]
%       lockedandcreep
%
%       For e.g. Topo=3 the Williams&Wadge topographic correction is carried out with the 
%       source depth adjusted according to the elevation (in 3 intervals in this case)
%
%  FA, June 2007   
%  FA displacement2rangechange included 07/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Put in the fixed parameters:

% direct assignments is faster than with a for loop
    N_disloc         = modelopt.N_disloc;
    N_fault          = modelopt.N_fault;
    N_mogi           = modelopt.N_mogi;
    N_penny          = modelopt.N_penny;
    N_mctigue        = modelopt.N_mctigue;
    N_yang           = modelopt.N_yang;
    N_squaredisloc   = modelopt.N_squaredisloc;
    N_multidisloc    = modelopt.N_multidisloc;
    N_visco1d        = modelopt.N_visco1d;
    N_lockedandcreep = modelopt.N_lockedandcreep ;
    N_peas           = modelopt.N_peas ;
    Topo             = modelopt.Topo;
    nu               = modelopt.dislocopt.nu    ;
    Layers           = modelopt.Layers;
   
    if Layers >0    YoungMod_layers = modelopt.layeropt.YoungMod_layers; end
    if Layers >0    mu_layers = YoungMod_layers./(2*(1+nu)); end
    if Layers >0    lam_layers = 2*mu_layers.*nu/(1-2*nu); end
    if Layers >0    scaleN           = modelopt.layeropt.scaleN; end
    
    
    u=[]; 
    ut=zeros(3,size(coord,2));
    if N_disloc >= 1   if Topo  ut=disloc1topo(par(1:10,:),coord,nu,hgt,Topo); elseif Layers ut=dislocLayers(par(1:10,:),coord,Layers,mu_layers,lam_layers,scaleN); else ut=disloc1(par(1:10,:),coord,nu); end; u=[u ut(:)]; par(1:10,:)=[]; end
    if N_disloc >= 2   if Topo  ut=disloc1topo(par(1:10,:),coord,nu,hgt,Topo); elseif Layers ut=dislocLayers(par(1:10,:),coord,Layers,mu_layers,lam_layers,scaleN); else ut=disloc1(par(1:10,:),coord,nu); end; u=[u ut(:)]; par(1:10,:)=[]; end
    if N_disloc >= 3   if Topo  ut=disloc1topo(par(1:10,:),coord,nu,hgt,Topo); elseif Layers ut=dislocLayers(par(1:10,:),coord,Layers,mu_layers,lam_layers,scaleN); else ut=disloc1(par(1:10,:),coord,nu); end; u=[u ut(:)]; par(1:10,:)=[]; end
    if N_disloc >= 4   if Topo  ut=disloc1topo(par(1:10,:),coord,nu,hgt,Topo); elseif Layers ut=dislocLayers(par(1:10,:),coord,Layers,mu_layers,lam_layers,scaleN); else ut=disloc1(par(1:10,:),coord,nu); end; u=[u ut(:)]; par(1:10,:)=[]; end
    if N_disloc >= 5  
        for i=5:N_disloc
            if ~Topo ut=disloc1(par(1:10,:),coord,nu);   else ut=disloc1topo(par(1:10,:),coord,nu,hgt,Topo);   end; u=[u ut(:)]; par(1:10,:)=[]; end
    end
        lamba=3.22e10; mu=3.43e10;
    if N_fault  >= 1   if ~Topo ut=fault(par(1:9,:),coord,lamba,mu);   else disp('noidea');   end; u=[u ut(:)]; par(1:9,:)=[]; end
    if N_fault  >= 2   if ~Topo ut=fault(par(1:9,:),coord,lamba,mu);   else disp('noidea');   end; u=[u ut(:)]; par(1:9,:)=[]; end
    if N_fault  >= 3   if ~Topo ut=fault(par(1:9,:),coord,lamba,mu);   else disp('noidea');   end; u=[u ut(:)]; par(1:9,:)=[]; end
    if N_fault  >= 4   if ~Topo ut=fault(par(1:9,:),coord,lamba,mu);   else disp('noidea');   end; u=[u ut(:)]; par(1:9,:)=[]; end
    

    if N_mogi   >= 1   if Topo  ut=mogitopo(par(1:4),coord,nu,hgt,Topo); elseif Layers  ut=MogiLayers(par(1:4),coord,Layers,mu_layers,lam_layers,scaleN); else  ut=mogi(par(1:4),coord,nu); end; u=[u ut(:)]; par(1:4) =[]; end
    if N_mogi   >= 2   if Topo  ut=mogitopo(par(1:4),coord,nu,hgt,Topo); elseif Layers  ut=MogiLayers(par(1:4),coord,Layers,mu_layers,lam_layers,scaleN); else  ut=mogi(par(1:4),coord,nu); end; u=[u ut(:)]; par(1:4) =[]; end
    if N_mogi   >= 3   if Topo  ut=mogitopo(par(1:4),coord,nu,hgt,Topo); elseif Layers  ut=MogiLayers(par(1:4),coord,Layers,mu_layers,lam_layers,scaleN); else  ut=mogi(par(1:4),coord,nu); end; u=[u ut(:)]; par(1:4) =[]; end
    

    if N_penny  >= 1   if ~Topo ut=penny_source(par(1:5),coord);  else ut=penny_sourcetopo(par(1:5),coord,hgt,Topo);  end; u=[u ut(:)]; par(1:5) =[]; end
    if N_mctigue>= 1   if ~Topo ut=mctigue(par(1:5),coord,1,nu);  else ut=mctiguetopo(par(1:5),coord,1,nu,hgt,Topo);  end; u=[u ut(:)]; par(1:5) =[]; end
    if N_yang   >= 1   if ~Topo ut=yang_source(par(1:8),coord,nu);else ut=yang_sourcetopo(par(1:8),coord,nu,hgt,Topo);end; u=[u ut(:)]; par(1:5) =[]; end

    if N_squaredisloc   >= 1  
        ut=squaredisloc(par(1:9,:),coord,nu);  u=[u ut(:)]; par(1:9,:)=[]; 
    end
    if N_multidisloc   >= 1  
                  len=10+length(modelopt.multidislocopt.ind);
                  ut=multidisloc(par(1:len),coord,modelopt.multidislocopt);u=[u ut(:)]; par(1:len)=[]; 
    end

    if N_visco1d        >= 1  ut=visco1d(par(1:10),coord,modelopt.visco1dopt);                u=[u ut(:)]; par(1:10) =[] ; end
    if N_lockedandcreep >= 1  ut(1,:)=lockedandcreep(par(1:5),coord(1,:));                    u=[u ut(:)]; par(1:5)  =[] ; end
     
    if N_peas           >= 1  ut = peas(par,coord) ;                                     u=[u ut(:)]; par(1:8)  =[] ; end

%
    u_radarlook = displacement2rangechange(u,radarlook);         %creates displacement in radar LOS direction (in direction of radarlook vector, can be cartesian for GPS enu)
%

    if nargout==3
        u_sum=sum(u,2) ; 
    end
