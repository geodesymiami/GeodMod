%function [V,Vs]=intgrV(fi,psi,h,Wt,t)
function [V]=intgrV(fi,h,Wt,t)
% V,Vs - volume of crack, volume of surface uplift
% fi,psi: basis functions
% t: interval of integration
%large=1e10;
V = sum(Wt.*fi.*t);
%Vs = sum(Wt.*fi.*(t-h*(h-t)./(h^2+t.^2)));
%V1 = sum(Wt.*(fi.*Q(0,t,0,41)));
%V2 = sum(Wt.*(fi.*Q(0,t,large,41)));
%V=V2-V1;
%I1=sqrt(2)*h*t.*Qr(h,t,1)
%I2=(h*Qr(h,t,3)-t.*Qr(h,t,2))/sqrt(2)
%I3=(h*Qr(h,t,2)+t.*Qr(h,t,3))/sqrt(2)
%Vs = sum(Wt.*(fi.*(I1+h*I2)+psi.*(I1./t-I3)))
%Vs2 = sum(Wt.*(fi.*(Q(h,t,large,41)+h*Q(h,t,large,51))+...,
%              psi.*(Q(h,t,large,41)./t-Q(h,t,large,61))));
%Vs=Vs2-Vs1;
