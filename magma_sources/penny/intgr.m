function [Uz,Ur]=intgr(r,fi,psi,h,Wt,t);
% function [Uz,Ur]=intgr(r,fi,psi,h,Wt,t)
% Uz(r),Ur(r) - vertical and radial displacements
% fi,psi: basis functions
% t: interval of integration

[s1,s2]=size(r);
Uz=zeros(size(r));
Ur=Uz;

for i=1:s1
 for j=1:s2
  Uz(i,j)= sum(Wt.*(fi.*(Q(h,t,r(i,j),1)+h*Q(h,t,r(i,j),2))+...,
            psi.*(Q(h,t,r(i,j),1)./t-Q(h,t,r(i,j),3))));
  Ur(i,j)=sum(Wt.*(psi.*((Q(h,t,r(i,j),4)-h*Q(h,t,r(i,j),5))./t-...,
           Q(h,t,r(i,j),6)+h*Q(h,t,r(i,j),7))-...,
           h*fi.*Q(h,t,r(i,j),8)));
 end
end
