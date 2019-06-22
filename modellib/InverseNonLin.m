function [modelopt]=InverseNonLin(dataset,inverseopt,modelopt)
% InverseModel-  non-linear inversion       
%
% usage: [mhat,mstar]=InverseNonLin(dataset,inverseopt,objfuncopt)
%
%Input:   dataset       structure containing input data              
%         inverseopt    structure with options for DistOp
%         objfuncopt    optiond for objective function                     
%
%Output:  mhat          dislocations with opening
%         mstar         Parameter of planes              
%
% inverseopt   -  options for inversion routine:
% 'algorithm'    Inversion algorithm (anneal,gibbs,gradient)       [default 'anneal']
%                algorithm is followed by gradient unless  FollowGradient='off'
% 'rstate'       random generator initialization                   [default clock]
%
% 'annealopt'    parameter to set annealopt                        [default 0]
%                    0  for   [1  1 1 2.5 0 0 1]
%                    1  for   [2  2 1 2.5 0 0 1]
%                    2  for   [4  3 2 2.5 0 0 1]
%                    3  for   [4  3 4 2.5 0 0 1]
%                    4  for   [10 5 4 2   0 0 1]
%                    5  for   [25 7 6 2   0 0 1]
%                    for an explanation of the options : help anneal
%
% 'gibbsopt'           parameter for Gibbs Sampling                        
% 'gibbsopt.Tschedopt' parameter for cooling schedule 
%           N  or  [TcStart TcEnd TcRobustNum CoolNum CoolSweeps GsSweeps]  (default 1)
%
%           TcStart       start temperature for Tc estimation  (e.g. -3.5,-4.7)
%           TcEnd         end   temperature  (TcStart=TcEnd to set Tc for Gibbs Sampling ) 
%           TcRobustNum   Number of temperature steps for Tc estimation   
%           CoolNum       Number of temperature steps to cool to Tc for Gibbs Sampling 
%           CoolSweeps    Sweeps at each temperature for cooling to Tc  
%           GsSweeps      Sweeps at Tc (actual Gibbs Sampling)   
%
%           if N is given  the following default parameters are set
%           0   for    [-3.9  -3.9  1  1  10  10   ]    (very fast, for testing)      
%           1   for    [-3.5  -4.9  3  1  10  10   ]    (very rapid Tc estimation)
%           2   for    [-3.5  -4.9  5  2  100 1000 ]    (very rapid Tc estimation)
%
%           example N=2: first the critical temperature will be estimated
%           i.e. the cooling schedule will have a lendbstop in ForwardModel_forBasemap at 32gth of 500 (100 sweeps at each of the
%           5 temperatures between -3.5 and -4.9 (the number of objective function evaluations 
%           is 500*igrid=4000). After having selected Tc the actual cooling schedule prescribes 
%           1200 sweeps, 200 of which for cooling and 1000 for the actual Gibbs Sampling.
%           Cooling will be conducted in 2 steps starting at T=Tc+1.5 (1.5 is hardwired in gibbs.m) 
%           each with 100 sweeps. 
%
% 'gibbsopt.igrid'       grid density  at each sweep                    [default 8]
% 'gibbsopt.plot_bins'   binning parameter for plotting ppds            [default 10]
% 'gibbsopt.plot_ninterp'interpolation paramater for plotting 2D ppds   [default 100]
% 'gibbsopt.nsave'       number of sweeps for each save                 [default 8]
% 'gibbsopt.runs '       number of gibbs sampling runs                  [default 8]
%
% V1.0  Falk Amelung, September 2006
% V1.1  Falk Amelung, June 2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global basedir
defaultopt=struct(                                            ...
        'DoIt'              ,   'on'                    ,     ...
        'objfunc'           ,   'GenericObjectiveFunction',   ...
        'algorithm'         ,   'anneal'                ,     ...
        'annealopt'         ,    0                      ,     ...
        'gridsearchopt'     ,    3                      ,     ...
        'FollowGradient'    ,   'off'                   ,     ...
        'fname'             ,   'model'                 ,     ...
        'PhaseRamp'         ,   'Ramp'                  ,     ...
        'FactorLin'         ,   'off'                   ,     ...
        'FactorNonLin'      ,   'off'                   ,     ...
        'FactorNonLinDelta' ,    0.2                    ,     ...
        'startmodel'        ,   'off'                   ,     ...
        'rstate'            ,    sum(100*clock)         ,     ...
        'QuickStop'         ,   'off'                   ,     ...
        'Plot'              ,   'off'      )            ;
defaultopt.gibbsopt=struct(                                   ... 
        'Tschedopt'         ,   'off'                   ,     ...
        'igrid'             ,    8                      ,     ...
        'matrix'            ,    0                      ,     ...
        'plot_bins'         ,    20                     ,     ...
        'plot_ninterp'      ,    100                    ,     ...
        'nsave'             ,    1000                   ,     ...
        'PlotLinearPar'     ,    'on'                   ,     ...
        'robustTc'          ,    0                      ,     ...
        'runs'              ,    1                      )     ;

global dir_out
[inverseopt]=process_defaultoptions(inverseopt,defaultopt);  %display(inverseopt)
f=fieldnames(inverseopt) ; for i=1:length(f) eval([char(f{i}) '= inverseopt.(f{i}) ;' ]) ; end

out_name=fullfile(dir_out,fname); 
if  ~DoIt                                    return; end
if  CheckInOut('',out_name)  
    load(out_name); 
    return; 
end     % return if model.mat exist
if  exist('modelopt.par','var')              return; end     % return if model parameters given
%if  modelopt.par_xy(1) || modelopt.par_lola(1) return; end     % return if model parameters given

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rand('state',rstate);
inverseopt.out_name=out_name ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Make bounds, identify free and fixed parameters, modify bounds, parnames according to SARmul, etc  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[inverseopt]   = InitializeInverseopt(dataset,inverseopt,modelopt);
[dataset]      = InitializeDataset(dataset,inverseopt) ;             %adds G_phaseramp to dataset
 objfunc       = inverseopt.objfunc    ;
 objfuncopt    = inverseopt.objfuncopt ;
 bounds        = inverseopt.bounds;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  ANNEAL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch algorithm
case{'Anneal'}
  if annealopt==0  annealopt=[1  1 1 2.5 0 0 1];  end
  if annealopt==1  annealopt=[2  2 1 2.5 0 0 1];  end
  if annealopt==2  annealopt=[4  3 2 2.5 0 0 1];  end
  if annealopt==3  annealopt=[4  3 4 2.5 0 0 1];  end
  if annealopt==4  annealopt=[10 5 4 2   0 0 1];  end
  if annealopt==5  annealopt=[25 7 6 2   0 0 1];  end
  logmessage(['annealopt: ' num2str(annealopt,' %g ')])

  sqrflag=true;
  estimate_inversion_time(objfunc,bounds,'Anneal',annealopt,dataset,objfuncopt,sqrflag) ;
  tic
  [mhat,f,models,energy,count]=anneal(objfunc,bounds,annealopt,dataset,objfuncopt,sqrflag);
  t=toc ; logmessage(sprintf('elapsed time for simulated annealing %.1f min',t/60))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   GIBBS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'Gibbs'}
      sqrflag=true;
  %
  % set gibbsopt according to input and initialize save file
  %
      f=fieldnames(gibbsopt) ; for i=1:length(f) eval([char(f{i}) '= gibbsopt.(f{i}) ;' ]) ; end
                            s   = 'TcStart TcEnd TcRobustNum CoolNum CoolSweeps GsSweeps'  ;
      if Tschedopt==0  Tschedopt=[ -3.9    -3.9      1          1      10       10   ];   end
      if Tschedopt==1  Tschedopt=[ -3.5    -4.9      3          1      10       10   ];   end
      if Tschedopt==2  Tschedopt=[ -3.5    -4.9      5          2     100     1000   ];   end
      if Tschedopt==3  Tschedopt=[ -3.5    -4.9      6          5     100     2000   ];   end
      if Tschedopt==4  Tschedopt=[ -3.5    -4.9      6         10     100    10000   ];   end
      if Tschedopt==5  Tschedopt=[ -3.5    -4.9      6         10     100   100000   ];   end
      if length(Tschedopt)~=6  errordgl('Wrong length of Tschedopt' - exiting'); end
         s=str2list(s); 
      for i=1:length(Tschedopt) defaultgibbsopt.(s{i})=Tschedopt(i); end 
      [gibbsopt]=process_defaultoptions(gibbsopt,defaultgibbsopt);  

      %if PlotLinearPar
      %    ind               = sort([inverseopt.freeind; inverseopt.linearind]);
      %else
      %    ind               = [inverseopt.freeind];
      %end
      
      ind                   = [inverseopt.freeind];
      
      gibbs_complete_flag   = false;
      gibbsopt.ParNames     = inverseopt.ParNames(ind);
      gibbsopt.ParForm      = inverseopt.ParForm (ind); 
      gibbsopt.ParNamesForm = inverseopt.ParNamesForm (ind); 
      if ~rstate  gibbsopt.rstate=sum(100*clock);   end
      gibbsopt.sfile        = [out_name '_gibbs'] ;

      if exist([gibbsopt.sfile '.mat'],'file')      % remove models_gibbs if too few models because of trying
          load(gibbsopt.sfile,'models','gibbsopt');
          if length(models)<length(find(gibbsopt.Tsched~=gibbsopt.Tsched(end)))
            error('#### Intentional interruption: models_gibbs.mat too small, consder gclean gibbs ######');   % gibbs sampling file may not have enough models
          end
      end

      if ~exist([gibbsopt.sfile '.mat'],'file')      % gibbs sampling to generate model_gibbs if solution file does not exist

          save(gibbsopt.sfile,'inverseopt') ;        %initialize sfile to append models,etc in gibbs (inverseopts not available in gibbs, bounds needed for plotting)

          logmessage(['Tschedopt: ' num2str(Tschedopt,' %g ')])

      %
      % compute Tc using Basu&Frazer method unless TcStart=TcEnd (i.e. Tc given)
      %
          estimate_inversion_time(objfunc,bounds,'Gibbs',gibbsopt,dataset,objfuncopt,sqrflag)

          tic
          [gibbsopt.Tc]          = calc_Tc(gibbsopt,objfunc,bounds,dataset,objfuncopt,sqrflag);
          t=toc ; logmessage(sprintf('elapsed time for Tc estimation %.1f min',t/60))
          [gibbsopt.Tsched,junk] = generate_Tsched(gibbsopt) ;

      %
      % Gibbs Sampling
      %
          tic
          [mhat,models,energy]=gibbs(objfunc,bounds,gibbsopt,dataset,objfuncopt,sqrflag);
          t=toc ; logmessage(sprintf('elapsed time for Gibbs Sampling %.1f min',t/60))
          gibbs_complete_flag = true;
      else
          load(gibbsopt.sfile); 
      end
           
      [gibbs_mean, gibbs_sigma] = normfit(models'); %fits gaussian to ppd's
       modelopt.par.gibbs_mean  = gibbs_mean';
       modelopt.par.gibbs_sigma = gibbs_sigma';
      
       if ~gibbs_complete_flag
          modelopt.par.gibbs_mean  = modelpar2invpar(modelopt.par.gibbs_mean, objfuncopt,-1);
          modelopt.par.gibbs_sigma = modelpar2invpar(modelopt.par.gibbs_sigma,objfuncopt,-1);
      end
      inverseopt.objfuncopt.modelopt                     = modelopt;      

  % Plot Probability distributions
  %
      gibbsopt.WHAT='PPD1D';
      out_name1=[gibbsopt.sfile '_ppd1d'];
      logplot('plot_gibbs',out_name1,models,energy,bounds,gibbsopt,inverseopt,objfuncopt);

      gibbsopt.WHAT='PPD2D';
      out_name2=[gibbsopt.sfile '_ppd2d'];
      logplot('plot_gibbs',out_name2,models,energy,bounds,gibbsopt,inverseopt,objfuncopt);
      
      if ~exist('mhat','var') error('#### Intentional interruption. Gibbs sampling underway. Dont worry! ######'); end      % only for geodmod runs for plotting of ppds  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  GridSearch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case{'GridSearch'}
     logmessage(['gridsearchopt: ' num2str(gridsearchopt,' %g ')])

     sqrflag=true;
     estimate_inversion_time(objfunc,bounds,'GridSearch',gridsearchopt,dataset,objfuncopt,sqrflag) ;
     tic
     sqrflag=1;
     [mhat,energy,models]=gridsearch(objfunc,bounds,gridsearchopt,dataset,objfuncopt,sqrflag);
     t=toc ; logmessage(sprintf('elapsed time for gridsearch %.1f min',t/60))

   % Plot Probability distributions

     gridsearchplotopt.WHAT='PPD2D';
     out_name2=['gridsearch_ppd2d'];
     gridsearchplotopt.plot_bins=10;
     gridsearchplotopt.modelopt=modelopt;
     gridsearchplotopt.ParNames     = inverseopt.ParNames(inverseopt.freeind); 
     gridsearchplotopt.ParForm      = inverseopt.ParForm (inverseopt.freeind); 
     gridsearchplotopt.ParNamesForm = inverseopt.ParNamesForm(inverseopt.freeind); 
     gridsearchplotopt.Tsched=1;    % needed so that plot_gibbs does not fail
     gridsearchplotopt.plot_ninterp=20;    % needed so that plot_gibbs does not fail
     %logplot('plot_gridsearch',out_name2,models,energy,bounds,gridsearchplotopt,inverseopt,objfuncopt);
     % plot_gridsearch is currently a copy of plot_gibbs. Needs to be adjusted to plot energy instead of
     % of the ppds.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  GRADIENT: LSQNONLIN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~isempty(mhat)) &&  ( inverseopt.FollowGradient || strcmp('Gradient',algorithm))
   if ~strmatch(algorithm,{'Anneal' 'Gibbs' 'GridSearch'})       % select startmodel if no annealing or gibbs
      if inverseopt.startmodel  
         mhat=inverseopt.startmodel ; 
      else
        looseind=find(bounds(:,1)~=bounds(:,2));
        mhat=bounds(looseind,1) + (bounds(looseind,2)-bounds(looseind,1))/2 ;
      end
   end

   sqrflag=false;
   lsqnonlinopt = optimset('Display','iter','MaxFunEvals',10000,'TolX',1e-6); 
   [mstar,F,r,exitflag,output,lambda,jacobian]=lsqnonlin(objfunc,mhat,bounds(:,1),bounds(:,2),lsqnonlinopt,dataset,objfuncopt,sqrflag);
else
   mstar=mhat;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  Convert solution to LatLong, save as structure  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mstar_full_xy                  = modelpar2invpar(mstar,objfuncopt,-1);
[j,j,j,j,mstar_full_xy,j,j,mlin]= feval(objfunc,mstar,dataset,objfuncopt,sqrflag);       % fix to retrieve linear model parameters
mstar_full_lola                 = modelpar_lola2xy(mstar_full_xy,plotdataopt.basemap,modelopt,-1);
modelopt.par.xy                 = mstar_full_xy;
modelopt.par.lola               = mstar_full_lola;
modelopt.Unit                   = dataset(1).Unit ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  Plot and save  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(out_name,'modelopt')

inverseopt.plotdataopt.Surface3D = false;
inverseopt.plotdataopt.modelopt  = modelopt ;
inverseopt.objfuncopt.modelopt   = modelopt ;

PrepareDataModelResidual_local_AndPlot(dataset,modelopt,inverseopt);

if QuickStop;   disp('### QuickStop -- that''s all, folks ###'); return; end     % Skip part of the plotting

PrepareDataModelResidual_lola_AndPlot   (dataset,modelopt,inverseopt);
