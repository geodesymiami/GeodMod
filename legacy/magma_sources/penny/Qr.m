function [K]=Qr(h,t,n)
% Kernels calculation 
K=[];
h2=h^2;
t2=t.^2;
h4=h^4;
t4=t.^4;
h6=h^6;
t6=t.^6;
h8=h^8;
t8=t.^8;
ht2=h2+t2;

switch n
case 1	%Q1
 K=2*(h6-h2*t4+t4*h2+h4*(2*t2+h2))./((2*h2)^(1.5)*(h4+t4-t2.*ht2+h2*(3*t2+h2)));

case 2	%Q2
 K=8*(h8+t8+h6*ht2-h4*t2.*ht2+h2*t4.*ht2-t6.*ht2)./(ht2*(2*h2)^(3.5));

case 3	%Q3
 rad=sqrt(4*h2*t2+(h2-t2).^2);
 K=2*(-h6-h4*t2+h2*t4+t6+(h4+t4).*rad)./(ht2.^2.*(t2-h2+rad).^(1.5));

end %switch n
