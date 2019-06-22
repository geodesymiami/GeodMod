function [str]=MakeStringForPlot(modelopt,inverseopt)
%MakeStringForPlot   - generates a string with solution, inverse modeling parameters 
%
%usage: [str_short,str_long]=MakeStringForPlot(par,inverseopt);
%
%Input:   same input list as PlotModel (for simplicity)  
%         rms_str    string tio be appended (e.g. with rms_values)
%
%Output:  plotstr    string with solution, etc  
%
%FA, September, 2005 

par        = modelopt.par.xy;
objfuncopt = inverseopt.objfuncopt;
N_disloc   = floor(length(par)/10);
bounds     = modelpar2invpar(inverseopt.bounds,objfuncopt,-1)    ;
if length(par)~=length(bounds) par=modelpar2invpar(par,objfuncopt,-1) ; end

rigidity = 3.0e10;              % Definition for seismic moment estimation

s=[] ;
if isfield(objfuncopt,'objfunc')        s=[s                      num2str(objfuncopt.objfunc)                  ] ; end
%TODO: New string depending on N_disloc,N_mogi,...
if objfuncopt.modelopt.N_disloc         s=[s ', dislocations: '   num2str(objfuncopt.modelopt.N_disloc)        ]; end
if objfuncopt.modelopt.N_squaredisloc   s=[s ', squaredisloc:'    num2str(objfuncopt.modelopt.N_squaredisloc)  ]; end
if objfuncopt.modelopt.N_multidisloc    s=[s ', multidisloc:'     num2str(objfuncopt.modelopt.N_multidisloc)   ]; end
if objfuncopt.modelopt.N_mogi           s=[s ', mogi: '           num2str(objfuncopt.modelopt.N_mogi)          ]; end
if objfuncopt.modelopt.N_mctigue        s=[s ', mcTigue: '        num2str(objfuncopt.modelopt.N_mctigue)       ]; end
if objfuncopt.modelopt.N_pCDM           s=[s ', pCDM: '           num2str(objfuncopt.modelopt.N_pCDM)          ]; end
if objfuncopt.modelopt.N_penny          s=[s ', penny: '          num2str(objfuncopt.modelopt.N_penny)         ]; end
if objfuncopt.modelopt.N_yang           s=[s ', yang: '           num2str(objfuncopt.modelopt.N_yang)          ]; end
if objfuncopt.modelopt.N_lockedandcreep s=[s ', lockedandcreep: ' num2str(objfuncopt.modelopt.N_lockedandcreep)]; end
if objfuncopt.modelopt.N_visco1d        s=[s ', visco1d: '        num2str(objfuncopt.modelopt.N_visco1d)       ]; end
if objfuncopt.modelopt.Topo             s=[s ', Topo: '           num2str(objfuncopt.modelopt.Topo)            ]; end
if inverseopt.algorithm                 s=[s ', algorithm: '      num2str(inverseopt.algorithm)                ]; end
if inverseopt.FollowGradient            s=[s ', FollowGradient: ' num2str(inverseopt.FollowGradient)           ]; end
if objfuncopt.FactorLin                 s=[s ', FactorLin: '      num2str(objfuncopt.FactorLin)                ]; end
if objfuncopt.FactorNonLin              s=[s ', FactorNonLin: '   num2str(objfuncopt.FactorNonLin)             ]; end
if objfuncopt.PhaseRamp                 s=[s ', PhaseRamp: '      num2str(objfuncopt.PhaseRamp)                ]; end
%if inverseopt.out_name                  s=[s ', saved as: '       num2str(inverseopt.out_name)                ]; end
str=s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind=1:length(par);
if length(par)>50 ind=1:50; end     % limit to display 50 model parameters

  strParNamesForm = list2str(objfuncopt.modelopt.ParNamesForm(ind)); 
  strParForm      = list2str(objfuncopt.modelopt.ParForm(ind));

  ctmpForm                        = objfuncopt.modelopt.ParForm(ind);
  ctmpForm(inverseopt.linearind)  = {objfuncopt.modelopt.ParNamesForm{inverseopt.linearind}};
  ctmpForm(inverseopt.fixind)     = {objfuncopt.modelopt.ParNamesForm{inverseopt.fixind}};
  tmpForm                         = list2str(ctmpForm);
  tmpbounds1                      = num2cell(bounds(ind,1));
  tmpbounds1(inverseopt.fixind)   = {inverseopt.ParType{inverseopt.fixind}};
  tmpbounds1(inverseopt.linearind)= {inverseopt.ParType{inverseopt.linearind}};

  str = strvcat( str, sprintf(strParNamesForm,objfuncopt.modelopt.ParNames{ind}) );
  str = strvcat( str, sprintf(tmpForm,        tmpbounds1{ind})                   );
  str = strvcat( str, sprintf(strParForm,     bounds(ind,2))                     );                     
  str = strvcat( str, sprintf(strParForm,     par(ind))                          );

  %if isfield(modelopt.par,'lola')
    %str = strvcat( str, sprintf(strParForm,     modelopt.par.lola(ind)) );
  %end

if isfield(modelopt.par,'gibbs_mean')
  str = strvcat( str, '  Gibbs sampling: mean and sigma for Gaussian:');
  str = strvcat( str, sprintf(strParForm,     modelopt.par.gibbs_mean(ind)) );
  str = strvcat( str, sprintf(strParForm,     modelopt.par.gibbs_sigma(ind)));
end

%s=sprintf(strParNamesForm,inverseopt.ParType{ind});          str2=char(str2,s) ;  % FA 7/08  taken out because of tmpnounds1
%form=[objfuncopt.modelopt.ParNamesForm{ind}]; s=sprintf(form,objfuncopt.modelopt.ParNames{ind})  %FA 7/08 works without list2str: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% if dislocation sources calculate geometric, seismic moments, Mw and opening volumes
%
if inverseopt.objfuncopt.modelopt.N_multidisloc >= 1
    x_unit                                  = 'km';                  % FA 2/2010. Don't know why to put this it is in multidisloc
    [disloc_par]                            = multidislocpar2dislocpar(par,inverseopt.objfuncopt.modelopt.multidislocopt,x_unit);
    inverseopt.objfuncopt.modelopt.N_disloc = length(disloc_par)/10;
    N_disloc                                = inverseopt.objfuncopt.modelopt.N_disloc;   %FA 2/2010: This is confusing. Should use N_disloc throughout ?
    par                                     = disloc_par;
end

if inverseopt.objfuncopt.modelopt.N_disloc >= 1
    
   len_ind             = 10*([0:N_disloc-1])+1 ;
   wid_ind             = 10*([0:N_disloc-1])+2 ;
   ss_ind              = 10*([0:N_disloc-1])+8 ;
   ds_ind              = 10*([0:N_disloc-1])+9 ;
   op_ind              = 10*([0:N_disloc-1])+10;
   slip_ind            = [ss_ind ds_ind op_ind];

   slipmag             = sqrt( par(ss_ind).^2 + par(ds_ind).^2 );
   geometric_moment    = sum(slipmag    .*par(len_ind).*par(wid_ind))*1e6;  % unit m^3 (*1e6 because len,wid is in km)
   ss_geometric_moment = sum(abs(par(ss_ind)).*par(len_ind).*par(wid_ind))*1e6;
   ds_geometric_moment = sum(abs(par(ds_ind)).*par(len_ind).*par(wid_ind))*1e6; 
   ss_percentage       = ss_geometric_moment/(ss_geometric_moment+ds_geometric_moment) * 100;
   ds_percentage       = ds_geometric_moment/(ss_geometric_moment+ds_geometric_moment) * 100;
   seismic_moment      = rigidity*geometric_moment
   
   op_volume           = sum(par(op_ind).*par(len_ind).*par(wid_ind))*1e6;   
   Mw                  = 2*log10(seismic_moment)/3 - 6.0;
   
   ss_max              = max(abs(par(ss_ind)));  %FA 2/2010 Don't know why the absolute was taken. Does not work well for Sjonni's Haiti solution
   ds_max              = max(abs(par(ds_ind)));
   op_max              = max(abs(par(op_ind)));
   %ss_max              = max(par(ss_ind));
   %ds_max              = max(par(ds_ind));      %FA 2/2010 If not the max is taken slipflag can be zero and no summary info is displayed. This seems to be a problem only for Sjonni's solution where slip is negarive
   %op_max              = max(par(op_ind));
   ss_min              = min(par(ss_ind));
   ds_min              = min(par(ds_ind));
   op_min              = min(par(op_ind));
   
   slip_flag           = logical(ss_max+ds_max);
   op_flag             = logical(op_max);
      
   if slip_flag str = strvcat( str, sprintf('Minimum displacement   %5.1f m strike-slip, %5.1f m disp-slip',  ss_min, ds_min             )); end
   if slip_flag str = strvcat( str, sprintf('Maximum displacement   %5.1f m strike-slip, %5.1f m disp-slip',  ss_max, ds_max             )); end             
   if op_flag   str = strvcat( str, sprintf('Maximum opening:       %5.1f m',                                 op_max                     )); end
   if slip_flag str = strvcat( str, sprintf('Fault type (moment %%): %5.1f %% strike-slip, %5.1f %% dip-slip',ss_percentage,ds_percentage)); end
   if slip_flag str = strvcat( str, sprintf('Moment magnitude:  Mw = %4.2f',                                  Mw                         )); end
   if op_flag   str = strvcat( str, sprintf('Opening volume:        %5.1f 10^6 m^3',                          op_volume*1e-6             )); end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

