function [FactorNonLin,sourceend]=GetFactorNonLin(par,dataset,objfuncopt);
%
% GetFactorNonLin   -  calculates an vector for multiplication with pred to account for factor differences between datasets      
%
%
% Input:   dataset   
%
% Output:  sourceend - index for last source parameter
%

FactorNonLin=objfuncopt.FactorNonLin ;

tmp=cell(length(dataset),1); [tmp{:}]=deal(dataset.Ndata); Ndata =cell2mat(tmp); datind=cumsum(Ndata);
tmp=cell(length(dataset),1); [tmp{:}]=deal(dataset.DataSet); DataSet=tmp;
tmp=cell(length(dataset),1); [tmp{:}]=deal(dataset.ind); ind=(tmp);

ieast   =strmatch('GPSeast' ,DataSet);
if isempty(ieast) lastSAR=length(DataSet)-1; else lastSAR=ieast-1; end              % lastSAR is actually length(DataSet), however -1 is necessary because there is not GPS data set

Factor=ones(datind(end),1);
if strcmp('SAR',FactorNonLin)
   fac=par(end);
   ind=[1:datind(lastSAR)]; Factor(ind)=Factor(ind)*fac ;
   %FactorNonLin=[fac*ones(datind(lastSAR),1); ones((datind(end)-datind(lastSAR)),1)] ;
   sourceend=length(par)-1;
elseif strcmp('GPS',FactorNonLin) 
   fac=par(end);
   ind=[datind(lastSAR)+1:length(Factor)]; Factor(ind)=Factor(ind)*fac ;
   sourceend=length(par)-1;
elseif strcmp('SARmul',FactorNonLin) 
   fac=par(end-lastSAR+1:end);
   i=1; if length(fac)>=i Factor(ind{i})=fac(i)*Factor(ind{i}); end
   i=2; if length(fac)>=i Factor(ind{i})=fac(i)*Factor(ind{i}); end
   i=3; if length(fac)>=i Factor(ind{i})=fac(i)*Factor(ind{i}); end
   i=4; if length(fac)>=i Factor(ind{i})=fac(i)*Factor(ind{i}); end
   i=5; if length(fac)>=i Factor(ind{i})=fac(i)*Factor(ind{i}); end
   sourceend=length(par)-length(fac);
end
FactorNonLin=Factor;

