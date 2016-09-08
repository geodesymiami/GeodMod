function [mhat,models,energy]=gibbs(FUN,bounds,gibbsopt,varargin)
%   gibbs sampling by Falk Amelung    - simulated annealing by Peter Cervelli
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
%'OPTIONS' specifies a number of annealing options (empty matrix or zeros for defaults):
%
%     gibbsoptTc.igrid   = grid spacing (default = 4).  Higher numbers permit finer levels
%                    of parameter discretization.
%     gibbsoptTc.matrix = flag that tells the algorithm whether the objective function
%                    can accept matrix input (default = 0); set to '1' if yes.  The
%                    objective function should accept models stored columnwise and
%                    return a vector of costs.  Writing the objective function this
%                    way can increase algorithm speed by 15-20%.
%     gibbsoptTc.Tsched = minimum T (default = 4)  
%     gibbsoptTc.maxT   = maximum T (default = 4) 
%     gibbsoptTc.numT   = number of Ts 
%     gibbsoptTc.numT   = number of sweeps 
%     gibbsoptTc.runs   = number of short runs 
%
%OUTPUTS:
%
%'model' is a matrix containing the bestmodels after each sweep.
%
%'energy' is a vector containing the costs corresponding to the models in 'model'.
%
%'T' are the temperatures.
%
% FA may 2002, based on code from Peter Cervelli.

%Check argument syntax

if nargin<2
	error('Usage: [mhat,F,model,energy,count]=anneal(FUN,bounds,OPTIONS,varargin)')
end

if nargin<3
	OPTIONS=[];
end

if size(bounds,2)~=2
	error('Second argument must be an nx2 matrix of parameter bounds, where n is the number of parameters.');
end

%Check OPTIONS

	if isempty(gibbsopt)
		error('gibbs:  needs gibbsopt');
    else
	    igrid        = gibbsopt.igrid  ;
		matrix       = gibbsopt.matrix;
		Tsched       = gibbsopt.Tsched;
		%rstate      = gibbsopt.rstate  ;
		robustTc     = gibbsopt.robustTc  ;
		runs         = gibbsopt.runs  ;
		nsave        = gibbsopt.nsave ;
		sfile        = gibbsopt.sfile ;
		ParNames     = gibbsopt.ParNames;
		ParForm      = gibbsopt.ParForm;
		ParNamesForm = gibbsopt.ParNamesForm;
        PlotLinearPar= gibbsopt.PlotLinearPar;
     end		

%Check bounds to make sure they're ok

	if max(bounds(:,1)>bounds(:,2))
		error('All the values in the first column of bounds must be less than those in the second.');
	end

  strNamesForm=list2str(ParNamesForm); strParForm=list2str(ParForm); 
  names=sprintf(strNamesForm,ParNames{:}); 

%Define constants

	p=size(bounds,1);
    delta=abs((bounds(:,1)-bounds(:,2)));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if (robustTc)    
       fprintf('\n'); 
       fprintf('#############################################################\n'); 
       fprintf('############### Determing critical temperature ##############\n'); 
	   nT=length(unique(Tsched));
       fprintf('######  %d runs, each %d sweeps at %d different Temperatures ######\n',runs,length(Tsched)/nT,nT)
       fprintf('######  Gridspacing: %d                                      ######\n',igrid)
   else
       nsweeps=length(find(Tsched==Tsched(end))); nbefore=length(Tsched)-nsweeps ;
       fprintf('\n'); 
       fprintf('######################################################################\n'); 
       fprintf('##################### Gibbs Sampling ###############################\n'); 
       fprintf('######  Sweeps for cooling,sampling: %d %d  igrid: %d  ####\n',nbefore,nsweeps,igrid)
       fprintf('######  Save every %d Sweeps in %12s  ##################\n',nsave,sfile)
   end	   
   energy=zeros(length(Tsched),runs);
   models=zeros(p,length(Tsched),runs);                                                   

count=0;
for r=1:runs
    %rand('state',(r-1)+rstate);
    startmodel=rand(p,1).*(bounds(:,2)-bounds(:,1)) + bounds(:,1) ;  %start each run with the same model
    bestmodel=startmodel ;  %start at each T with the same model
	fprintf( ' #################### RUN: %d ###############################\n ',r);
    %fprintf('#### Sweep parameters: %s \n',names); 
    fprintf('    sweep  T   energy  %s \n',names); 
	fprintf( ['            startmodel: ' strParForm ],startmodel(:)); fprintf(['\n']);

    for i=1:length(Tsched)
      if (robustTc & i>1 & Tsched(i)~=Tsched(i-1)) bestmodel=startmodel ; end  %start at each T with the same model if searching Tc
      temp=Tsched(i);
         for x=1:p                   % Visit each parameter
	      	if delta(x)
            %Evaluate objective function at each permissible value
			vals=[0:1/igrid:(1-1/igrid)] + rand(1,igrid)/igrid;
           	v=[bestmodel(x) bounds(x,1)+vals*delta(x)];
           	modelmatrix=bestmodel*ones(1,length(v));
           	modelmatrix(x,:)=v;
			NM=size(modelmatrix,2);
			if matrix
				O=feval(FUN,modelmatrix,varargin{:});
			else
				O=zeros(NM,1);
				for e=1:NM 
					O(e)=feval(FUN,modelmatrix(:,e),varargin{:});					
	                count=count+1;
          		end
			end

          %Form exponential probability distribution
           
           	[dist,nanflag]=MakePDF(temp,O);
			if nanflag~=0
				for u=1:length(nanflag)
				   disp(['Warning: the cost function generated a NaN for the following model:']);
				   disp(modelmatrix(nanflag(u)))
				end
            end

	    	%Sample from probability distribution
           	s=find(cumsum(dist)>=rand);
			s=s(1);
           	bestmodel(x,1)=modelmatrix(x,s);
		    end
		 end       % end x-loop
         energy(i,r)=O(s);
       % now assign best model (after sweep over all model parameters) to models. 

	     models(:,i,r)=bestmodel;

         if (robustTc)   
             if (i==length(Tsched) | (i>1 & Tsched(i)~=Tsched(i+1)) )   
	            fprintf( [' endmodel:  %6.3f %7.5f  ' strParForm '\n'],log10(temp),energy(i,r),bestmodel(:));
		     end
	     else
	         if(mod(i,nsave/5)==0 | i==length(Tsched)) 
			    fprintf( ['    %4d %6.3f %7.5f ' strParForm],i,log10(temp),energy(i,r),bestmodel(:)); fprintf(['\n']);
	            if(mod(i,nsave) | i==length(Tsched)==0)   
				   tmpmodels=models;tmpenergy=energy;tmpTsched=Tsched;   
				   dellist=find(energy==0); energy(dellist)=[];models(:,dellist)=[];Tsched(dellist)=[];               % shorten allocated arrays 
				   tdellist=find(Tsched~=Tsched(end)); Tsched(tdellist)=[];models(:,tdellist)=[];energy(tdellist)=[]; % remove data not at Tc
				   if ~isempty(models) save(sfile,'models','energy','Tsched','bounds','gibbsopt','-append');   end
				   models=tmpmodels;energy=tmpenergy;Tsched=tmpTsched;
				end
	         end
         end
    end	          % end of T-loop
end               % end of runs loop

mhat=models(:,find(energy==min(energy(:)),1)); 

% FA 10/2008: the linear model parameters need to calculated separately because the non-linear ones are selected individually after a sweep
% (saving them from each model call does not work as bestmodel(x,1)=modelmatrix(x,s) selects one model parameter at a time

    logmessage('running objective function to calculate linear model parameters...')
    tmpmodels=zeros(p+length(varargin{2}.linearind)+length(varargin{2}.fixind),length(Tsched),runs);     % FA Feb 2008: calculate linear model params
    for r=1:runs
        for i=1:length(models)
            [j,j,j,j,tmpmodels(:,i,r)]=feval(FUN,models(:,i,r),varargin{:});
        end
    end
    %tmpmodels(varargin{2}.fixind,:) = [];
    models                          = tmpmodels;
    save(sfile,'models','energy','Tsched','bounds','gibbsopt','-append');     %save including linear model parameters
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pdf,nanflag]=MakePDF(temp,v)
%Forms exponential probability distribution given a 
% temperature and vector of costs.
%Internal function for simulated annealing algorithm.

bad=find(isnan(v));
if isempty(bad)
	pdf=eprob(temp,v);
	nanflag=0;
else
	good=find(~isnan(v));
	w=v(good);
	pdf=eprob(temp,w);
	pdf(good)=pdf;
	pdf(bad)=0;
	nanflag=bad;
end

function [pdf]=eprob(temp,v)
%Scales cost vector and calculates exponential probability distribution.  
%The scaling is necessary to permit wide ranges in temperature.
%Internal function for simulated annealing algorithm.

toobig=708.3964185322641;
pdf=v/temp;
mpdf=max(pdf);

if mpdf>toobig
   scale=mpdf/toobig;
   pdf=exp(-pdf/scale);
   pdf=pdf/max(pdf);
   pdf=pdf.^scale;
else
   pdf=exp(-pdf);
   pdf=pdf/max(pdf);
end

pdf=pdf/sum(pdf);
