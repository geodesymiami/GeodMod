function [mhat,F,model,energy,count]=anneal(FUN,bounds,OPTIONS,varargin)
%   anneal         - simulated annealing by Peter Cervelli
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
%     OPTIONS(1) = scale of cooling schedule (default = 4).  Higher numbers
%                  produce more exhaustive searches.
%     OPTIONS(2) = number of individual annealing runs (default = 3).  Higher 
%                  numbers produce more exhaustive searches and reduce dependency
%                  on correctly guessing critical temperature.
%     OPTIONS(3) = grid spacing (default = 4).  Higher numbers permit finer levels
%                  of parameter discretization.
%     OPTIONS(4) = temperature scale. To tweak this parameter, inspect
%                  a graph of 'energy'.  Values higher than 3 or lower than 1 are
%                  seldom, if ever, warranted.  The default tries several different
%                  value and works well for most problems.
%     OPTIONS(5) = flag that tells the algorithm whether the objective function
%                  can accept matrix input (default = 0); set to '1' if yes.  The
%                  objective function should accept models stored columnwise and
%                  return a vector of costs.  Writing the objective function this
%                  way can increase algorithm speed by 15-20%.
%     OPTIONS(6) = flag that tells the algorithm whether to try the Nelder-Mead
%                  simplex method (MATLAB's 'fmin' function) to improve answer
%                  (default = 0); set to '1' if yes.  Setting this flag to '2'
%                  causes the algorithm to use the Optimization Toolbox constrained
%                  optimization function, 'constr'.
%     OPTIONS(7) = flag the tells the algorithm whether to display informative
%                  output (default = 0); set to '1' if yes.
%
%OUTPUTS:
%
%'mhat' is the best model found.
%
%'F' is the cost associated with the best model.
%
%'model' is a matrix containing the bestmodels after each sweep.
%
%'energy' is a vector containing the costs corresponding to the models in 'model'.
%
%'count' is the total number of models checked.
%
%Version 1.0  Peter Cervelli 4-26-98.

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

	if isempty(OPTIONS)
		scale=4;
		runs=3;
		grid=4;
		ts=linspace(1.5,2.5,runs);
		matrix=0;
		newton=1;
		talk=1;
	else
		OPTIONS(8)=0;
		if OPTIONS(1)
			scale=OPTIONS(1);
		else
			scale=5;
		end

		if OPTIONS(2)
			runs=OPTIONS(2);
		else
			runs=2;
		end

		if OPTIONS(3)
			grid=OPTIONS(3);	
		else
			grid=4;
		end

		if OPTIONS(4);
			ts=ones(runs,1)*OPTIONS(4);
		else
			ts=linspace(1.5,2.5,runs);
		end

		matrix=OPTIONS(5);
		newton=OPTIONS(6);
		talk=OPTIONS(7);
	end

%Check bounds to make sure they're ok

	if max(bounds(:,1)>bounds(:,2))
		error('All the values in the first column of bounds must be less than those in the second.');
	end

%Define constants

	p=size(bounds,1);
	
	count=zeros(runs,1);
	energy=Inf;
     vals=[2.^-(1:grid)];
     vals=[vals,0,-vals];
     delta=0.5*abs((bounds(:,1)-bounds(:,2)));

%Loop through runs

for k=1:runs
c=0;

	bestmodel=rand(p,100).*((bounds(:,2)-bounds(:,1))*ones(1,100))+bounds(:,1)*ones(1,100);
	
	if matrix
		O=feval(FUN,bestmodel,varargin{:});
	else
		O=zeros(100,1);
		for e=1:100 
			O(e)=feval(FUN,bestmodel(:,e),varargin{:});					
     	end
	end
     tc=log10(mean(O))-ts(k);
	[v,i]=min(O);
     bestmodel=bestmodel(:,i);

if talk
	fprintf('\n\nBeginning run #%02d. Critical temperature at %3.4f.\n',k,10^tc);
	fprintf('------------------------------------------------\n\n');
	fprintf('f-Calls\t\tTemperature\tMinimum f-Value\n')
	fprintf('------------------------------------------------\n');
end

%Create cooling schedule from critical temperature and scale

	x=scale*[1 2 4 6 10 6 4 2 1];
	t=sum(x);
	temp=logspace(tc+1,tc-1,9);
	T=zeros(t,1);
	C=1;
	for i=1:9
		for j=1:x(i)
			T(C)=temp(i);
			C=C+1;
		end
	end

%Begin annealing

for w=1:t
   
	temp=T(w);
   	c=c+1;
	
	if talk
		if c/10==floor(c/10)
			fprintf('%7d\t\t%8.4f\t\t%8.3f\n',count(k),temp,min(energy(1:c-1,k)));
		end
	end
	%Visit each parameter
   
  	 for x=1:p
      
		if delta(x)

          %Evaluate objective function at each permissible value
           
           	v=bestmodel(x)+vals*delta(x);
          	v=v(find((v<=bounds(x,2))&(v>=bounds(x,1))));
           	modelmatrix=bestmodel*ones(1,length(v));
           	modelmatrix(x,:)=v;
			NM=size(modelmatrix,2);
          	count(k)=count(k)+NM;
			if matrix
				O=feval(FUN,modelmatrix,varargin{:});
			else
				O=zeros(NM,1);
				for e=1:NM 
					O(e)=feval(FUN,modelmatrix(:,e),varargin{:});					
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
           	energy(c,k)=O(s);
           	model(:,c)=bestmodel;
   		
		end
	end
end


[F(k,1),i]=min(energy(:,k));
mhat(:,k)=model(:,i);

if newton
	if newton==1
		mstar=fmins(FUN,mhat(:,k),[],[],varargin{:});
		Ostar=feval(FUN,mstar,varargin{:});
		if Ostar<F(k,1)
			if mstar>bounds(:,1) & mstar<bounds(:,2)
				mhat(:,k)=mstar;
				F(k,1)=Ostar;
				if talk
					fprintf('\nSimplex method lowered cost and remained within constraints.\n\n')
				end
			else
				if talk
					fprintf('\nSimplex method lowered cost but failed to remain within constraints.\n\n')
				end
			end
		end
	elseif newton==2
		mstar=constr(FUN,mhat(:,k),[],bounds(:,1),bounds(:,2),[],varargin{:});
		Ostar=feval(FUN,mstar,varargin{:});
		if Ostar<F(k,1)
             mhat(:,k)=mstar;
             F(k,1)=Ostar;
             if talk
                fprintf('\nConstrained optimization lowered cost.\n\n')
             end
		end	

	end
end

end

[F,i]=min(F);
mhat=mhat(:,i);


function [pdf,nanflag]=MakePDF(temp,v)
%Forms exponential probability distribution given a temperature and vector of costs.
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
%Scales cost vector and calculates exponential probability distribution.  The scaling
%is necessary to permit wide ranges in temperature.
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

