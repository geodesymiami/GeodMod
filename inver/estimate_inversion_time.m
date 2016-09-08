function []=estimate_inversion_time(objfunc,bounds,algorithm,opt,varargin);
%estimate_inversion_time    - 
%ANNEAL     [mhat,F,model,energy,count]=anneal(FUN,bounds,OPTIONS,x1,x2...,xn)
%
%Simulated annealing algorithm that tries to find a minimum to the function 'FUN'.
%
%INPUTS:
%
%'FUN' specifies the objective function.  This function should accept an input model vector
%as the first input a scalar cost. Additional context-specific argumentscan be passed to the
%objective function by passing them to 'ANNEAL' after 'OPTIONS'.
%
%'bounds' specifies the upper and lower limits that each model parameter can take on.  This matrix
%must have as many rows as model parameters and two columns.
%
% FA may 2007

%Check argument syntax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [Anneal, GridSearch, Gibbs]=deal(false);
 if strcmp('Anneal',    algorithm)    Anneal    = true; end
 if strcmp('Gibbs',     algorithm)    Gibbs     = true; end
 if strcmp('GridSearch',algorithm)    GridSearch= true; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logmessage(sprintf('benchmarking..'));
t_thresh=0.5;
par=bounds(:,1)+abs((bounds(:,1)-bounds(:,2)))/2;

    tic;
        O=feval(objfunc,par,varargin{:}); 
    t=toc; 
    t_call=t;

if t < t_thresh
    tic;
       for i=1:10   O=feval(objfunc,par,varargin{:});  end
    t=toc; 
    t_call=t/10;
end

if t < t_thresh
    tic;
       for i=1:100   O=feval(objfunc,par,varargin{:});  end
    t=toc ;
    t_call=t/100;
end

logmessage(sprintf('Elapsed time for 1 call of objective function: %.3f sec',t_call));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Gibbs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Gibbs
     f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end
     gibbsopt=opt;

     strNamesForm=list2str(ParNamesForm); strParForm=list2str(ParForm); 
     names=sprintf(strNamesForm,ParNames{:}); 

     [Tsched,Tsched_robust] = generate_Tsched(gibbsopt) ;
     gibbsopt.robustTc = true ;
     gibbsopt.Tsched   = Tsched_robust;
     gibbsopt.runs     = 5 ;
  
     %tmp_gibbsopt=gibbsopt;
     %tmp_gibbsopt.runs=1;
     %tmp_gibbsopt.Tsched=tmp_gibbsopt.Tsched(1);
     %tic 
     %[mhat,modelsTc,energyTc]=gibbs(objfunc,bounds,tmp_gibbsopt,varargin{:}); 
     %t_sweep=toc
     %logmessage(sprintf('Elapsed time for 1 sweep: %.3f sec',t_sweep));

     icall_sweep    = length(bounds)*(igrid+1);
     isweep_robust  = 5*TcRobustNum*CoolSweeps ;
     isweep_gs      = runs*CoolNum*CoolSweeps+GsSweeps ;
     icall_robust   = isweep_robust*icall_sweep;
     icall_gs       = isweep_gs    *icall_sweep;
     logmessage(sprintf('Estimated time for robust Tc estimation  (%6d sweeps,%6d f-calls): %.1f min',isweep_robust,icall_robust,icall_robust*t_call/60));
     logmessage(sprintf('Estimated time for gibbs sampling (%6d sweeps,%6d f-calls): %.1f min',isweep_gs    ,icall_gs    ,icall_gs    *t_call/60));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Anneal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Anneal
     annealopt=opt;
     scale = annealopt(1);                 % scale (higher numbers produce more exhaustive searches)
	 runs  = annealopt(2);                 % # of runs
     grid  = annealopt(3);                 % grdi spacing
     ts    = ones(runs,1)*annealopt(4) ;   %temperature scale

     x=scale*[1 2 4 6 10 6 4 2 1] ; t=sum(x);
     vals=[2.^-(1:grid)]; vals=[vals,0,-vals];

     icall_sweep = length(bounds)*length(vals);
     isweep_run  = t ;
     logmessage(sprintf('Estimated time for empirical Tc estimation  (100 f-calls): %.1f min',100*t_call/60));
     logmessage(sprintf('Estimated time for simulated annealing  (%d f-calls): %.1f min',runs*isweep_run*icall_sweep, runs*isweep_run*icall_sweep*t_call/60));
     logmessage(        'actual time will be 10-20% less because some candidate gridpoints are out of bounds' );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% GridSearch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if GridSearch
     n_gridpoints=(opt+1)^size(bounds,1);
     logmessage(sprintf('Estimated time for grid search on %d grid points: %.1f min',n_gridpoints,(t_call*n_gridpoints)/60));
end
