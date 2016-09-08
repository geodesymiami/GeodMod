function [Tsched,Tsched_robust] = generate_Tsched(opt) ;
%generate_Tsched      - generate cooling schedules based on gibbsopt
%
%FA, May, 2007 
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end

Tsched_robust=logspace(TcStart, TcEnd, TcRobustNum);
Tsched_robust=reshape(  repmat(Tsched_robust,CoolSweeps,1) ,1,CoolSweeps*length(Tsched_robust) );

%Tsched = logspace(Tc+1.5,Tc,CoolNum) ; 
if ~exist('Tc','var') Tc=TcStart ; end
Tsched = logspace(Tc+1.5,Tc,CoolNum+1) ; Tsched(end)=[]; 
tmp    = repmat(Tsched,CoolSweeps,1); 
Tsched = [tmp(:)' ones(1,GsSweeps)*10^Tc];
str=sprintf(' TcTstart %f TcEnd %f TcRobustNum %d CoolSweeps %d CoolNum=%d, GsSweeps=%d',TcStart,TcEnd,TcRobustNum,CoolNum,CoolSweeps,GsSweeps);   
%logmessage(str);

