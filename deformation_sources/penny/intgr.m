function [Uz,Ur]=intgr(r,fi,psi,h,Wt,t)
% function [Uz,Ur]=intgr(r,fi,psi,h,Wt,t)
% Uz(r),Ur(r) - vertical and radial displacements
% fi,psi: basis functions
% t: interval of integration

[s1,s2]=size(r);

% To speed up, use matrix operation instead of loop, Yunjun, 2015-11-10
rr   = repmat(reshape(r,s1*s2,1),size(t));
tt   = repmat(t,  s1*s2,1);
Wt2  = repmat(Wt, s1*s2,1);
fi2  = repmat(fi, s1*s2,1);
psi2 = repmat(psi,s1*s2,1);

Qf = Q(h,tt,rr,1:8);
Uz = sum(Wt2.*fi2.*(Qf(:,:,1) + h*Qf(:,:,2) + psi2.*Qf(:,:,1)./tt - Qf(:,:,3)),2);
Ur = sum(Wt2.*(psi2.*((Qf(:,:,4) - h*Qf(:,:,5))./tt - Qf(:,:,6) + h*Qf(:,:,7)) - h*fi2.*Qf(:,:,8)),2);
Uz = reshape(Uz,s1,s2);
Ur = reshape(Ur,s1,s2);


% Uz=zeros(size(r));
% Ur=Uz;
% for i=1:s1
%     for j=1:s2
%         Uz(i,j)= sum(Wt.*(fi.*(Qy(h,t,r(i,j),1)+h*Qy(h,t,r(i,j),2))+...,
%             psi.*(Qy(h,t,r(i,j),1)./t-Qy(h,t,r(i,j),3))));
%         Ur(i,j)= sum(Wt.*(psi.*((Qy(h,t,r(i,j),4)-h*Qy(h,t,r(i,j),5))./t-...,
%             Qy(h,t,r(i,j),6)+h*Qy(h,t,r(i,j),7))-...,
%             h*fi.*Qy(h,t,r(i,j),8)));
%     end
% end
