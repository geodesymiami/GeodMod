function [Uz]=intgrUz(r,fi,psi,h,Wt,t)
% Uz(r) - vertical displacements
% fi(t),psi(t): basis functions
% t: interval of integration 
[s1,s2]=size(r);
Uz=zeros(size(r));

for i=1:s1
 for j=1:s2
  Uz(i,j)= sum(Wt.*(fi.*(Q(h,t,r(i,j),1)+h*Q(h,t,r(i,j),2))+...,
            psi.*(Q(h,t,r(i,j),1)./t-Q(h,t,r(i,j),3))));
 end
end
