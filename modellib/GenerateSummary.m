function [plotstr] = GenerateSummary(modelopt,dataset,inverseopt)
%   GenerateSummary       - Calculates RMS, Weights, and writes them into a string                  
%   FA 7/2008
%
  %logmessage(sprintf('[]=%s',mfilename));
  
  objfuncopt = inverseopt.objfuncopt ;
  objfunc    = inverseopt.objfuncopt.objfunc ; 

  modelpar = modelopt.par.xy;
  invpar   = modelpar2invpar(modelpar,inverseopt.objfuncopt,1) ;
     
  sqrflag=true;
  %whos invpar
  [resi,pred,u,rms,npar,weights,rms_unitsig,mlin,dataset] = feval(objfunc,invpar,dataset,objfuncopt,sqrflag);         %FA Feb 2008. Needed to make it work for ProjectProfile
  par=npar;

  if        strcmp(dataset(1).Unit,'m/yr')  Unit='mm/yr';     % We assume that Unit is m/yr
     elseif strcmp(dataset(1).Unit,'m')     Unit='mm'   ;    % FA Oct 2007. Not sure whetehr correct
     else                                   Unit='?';
  end

  str=sprintf('TOTAL RMS: %5.2f mm  (',sqrt(resi/length(pred))*1000); 
                                for i=1:length(dataset) s=sprintf('%s %5.2f%s, ',dataset(i).DataSet,rms(i)*1000,Unit)        ;str =[str s ];end;str=[str ')']; 
  str1='NUMBER OF DATA POINTS ';for i=1:length(dataset) s=sprintf('%s %-10d, '  ,dataset(i).DataSet,dataset(i).Ndata)        ;str1=[str1 s];end;str=strvcat(str,str1); 
  str1='RMS FOR UNIT SIGMAS   ';for i=1:length(dataset) s=sprintf('%s %5.2f%s, ',dataset(i).DataSet,rms_unitsig(i)*1000,Unit);str1=[str1 s];end;str=strvcat(str,str1); 
  str1='WEIGHTS:              ';for i=1:length(dataset) s=sprintf('%s %5.2f'  ,dataset(i).DataSet,weights(i),'%    ,')        ;str1=[str1 s];end;str=strvcat(str,str1);  
  
  [plotstr]  = strvcat( str,     MakeStringForPlot(modelopt,inverseopt)                      );
  [plotstr]  = strvcat( plotstr, MakeStringLinearInversionInfo(inverseopt,dataset,mlin) );
