function [Tc] = calc_Tc(gibbsopt,objfunc,bounds,dataset,objfuncopt,sqrflag) ;
%calc_Tc      - calculate critical temperature for Gibbs Sampling
%
%FA, May, 2007 
f=fieldnames(gibbsopt) ; for i=1:length(f) eval([char(f{i}) '= gibbsopt.(f{i}) ;' ]) ; end
  if TcStart==TcEnd
     Tc=gibbsopt.TcStart ;
     str=sprintf('ManualTc: %6.3f' ,Tc);
  else
    [junk,Tsched_robust] = generate_Tsched(gibbsopt) ;
     gibbsopt.robustTc = true ;
     gibbsopt.Tsched   = Tsched_robust;
     gibbsopt.runs     = 5 ;

     [mhat,modelsTc,energyTc]=gibbs(objfunc,bounds,gibbsopt,dataset,objfuncopt,sqrflag);
     [Tc]=find_robustTc(modelsTc,energyTc,gibbsopt);
     out_name=[gibbsopt.sfile '_robustTc'];
     logplot('plot_robustTc',out_name,modelsTc,energyTc,gibbsopt)
     str=sprintf('RobustTc: %6.3f' ,Tc);
  end

  logmessage(str);



