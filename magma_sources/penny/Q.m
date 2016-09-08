function [K]=Q(h,t,r,n)
% Kernels calculation 
K=[];
E=h^2+r^2-t.^2;
D=(E.^2+4*h^2*t.^2).^(0.5);
D3=D.^3;
%i=sqrt(-1);

switch n
case 1	%Q1
 K=sqrt(2)*h*t./(D.*sqrt(D+E));

case 2	%Q2
 K=1/sqrt(2)./D3.*(h*sqrt(D-E).*(2*E+D)-t.*sqrt(D+E).*(2*E-D));

case 3	%Q3
 K=1/sqrt(2)./D3.*(h*sqrt(D+E).*(2*E-D)+t.*sqrt(D-E).*(2*E+D));

case 4	%Q4
 K=t/r-sqrt(D-E)/r/sqrt(2);
% K=imag(sqrt(E+2*i*h.*t)/r);

case 5	%Q5
 K=-(h*sqrt(D-E)-t.*sqrt(D+E))./D/r/sqrt(2);

case 6	%Q6
 K=1/r-(h*sqrt(D+E)+t.*sqrt(D-E))./D/r/sqrt(2);
% K=1/r-real((h-i*t)/r./sqrt(E-2*i*h.*t));

case 7	%Q7
 K=r*sqrt(D+E).*(2*E-D)./D3/sqrt(2);
% K=real(r.*(E-2*i*h.*t).^(-3/2));

case 8	%Q8
 K=r*sqrt(D-E).*(2*E+D)./D3/sqrt(2);
% K=imag(r.*(E-2*i*h.*t).^(-3/2));

case 41	%Q4*r
 K=t-sqrt(D-E)/sqrt(2);
% K=imag(sqrt(E+2*i*h.*t)/r);

case 51	%Q5*r
 K=-(h*sqrt(D-E)-t.*sqrt(D+E))./D/sqrt(2);

case 61	%Q6*r
 K=1-(h*sqrt(D+E)+t.*sqrt(D-E))./D/sqrt(2);

end %switch n
