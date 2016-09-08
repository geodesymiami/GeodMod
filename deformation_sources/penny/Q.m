function [K]=Q(h,t,r,n)
% Kernels calculation
% Support matrix 't' and array 'n' input, Yunjun, 2015-11-10

K=[];
E=h^2+r.^2-t.^2;
D=(E.^2+4*h^2*t.^2).^(0.5);
D3=D.^3;
%i=sqrt(-1);

% switch n
%     case 1	%Q1
%         K=sqrt(2)*h.*t./(D.*sqrt(D+E));
%         
%     case 2	%Q2
%         D3=D.^3;
%         K=1/sqrt(2)./D3.*(h.*sqrt(D-E).*(2.*E+D)-t.*sqrt(D+E).*(2.*E-D));
%         
%     case 3	%Q3
%         D3=D.^3;
%         K=1/sqrt(2)./D3.*(h.*sqrt(D+E).*(2.*E-D)+t.*sqrt(D-E).*(2.*E+D));
%         
%     case 4	%Q4
%         K=t./r-sqrt(D-E)./r/sqrt(2);
%         % K=imag(sqrt(E+2*i*h.*t)/r);
%         
%     case 5	%Q5
%         K=   -(h.*sqrt(D-E)-t.*sqrt(D+E))./D./r/sqrt(2);
%         
%     case 6	%Q6
%         K=1./r-(h.*sqrt(D+E)+t.*sqrt(D-E))./D./r/sqrt(2);
%         % K=1/r-real((h-i*t)/r./sqrt(E-2*i*h.*t));
%         
%     case 7	%Q7
%         D3=D.^3;
%         K=r.*sqrt(D+E).*(2.*E-D)./D3/sqrt(2);
%         % K=real(r.*(E-2*i*h.*t).^(-3/2));
%         
%     case 8	%Q8
%         D3=D.^3;
%         K=r.*sqrt(D-E).*(2.*E+D)./D3/sqrt(2);
%         % K=imag(r.*(E-2*i*h.*t).^(-3/2));
%         
%     case 41	%Q4*r
%         K=t-sqrt(D-E)/sqrt(2);
%         % K=imag(sqrt(E+2*i*h.*t)/r);
%         
%     case 51	%Q5*r
%         K= -(h.*sqrt(D-E)-t.*sqrt(D+E))./D/sqrt(2);
%         
%     case 61	%Q6*r
%         K=1-(h.*sqrt(D+E)+t.*sqrt(D-E))./D/sqrt(2);
%         
% end %switch n

%Support array 'n' input, Yunjun, 2015-11-10
sq_DE1 = sqrt(D+E);
sq_DE2 = sqrt(D-E);
[s1,s2] = size(t);  s3 = length(n);
if s3 > 1
    K=zeros(s1,s2,s3);
    for i=1:s3
        K(:,:,i) = Q_Kn(h,t,r,n(i),E,D,D3,sq_DE1,sq_DE2);
    end
elseif s3 == 1    % same output as Q if len(n)==1
    K = Q_Kn(h,t,r,n,E,D,D3);
else
    disp('Empty array n!')
    return
end

end

%% Support array 'n' input, Yunjun, 2015-11-10
function [Km] = Q_Kn(h,t,r,m,E,D,D3,sq_DE1,sq_DE2)

if m == 1;    Km=sqrt(2)*h.*t./(D.*sq_DE1);                                      end

if m == 2;    Km=1/sqrt(2)./D3.*(h.*sq_DE2.*(2.*E+D)-t.*sq_DE1.*(2.*E-D));    end

if m == 3;    Km=1/sqrt(2)./D3.*(h.*sq_DE1.*(2.*E-D)+t.*sq_DE2.*(2.*E+D));    end

if m == 4;    Km=t./r-sq_DE2./r/sqrt(2);                                         end

if m == 5;    Km=    -(h.*sq_DE2-t.*sq_DE1)./D./r/sqrt(2);                     end

if m == 6;    Km=1./r-(h.*sq_DE1+t.*sq_DE2)./D./r/sqrt(2);                    end

if m == 7;    Km=r.*sq_DE1.*(2.*E-D)./D3/sqrt(2);                                end

if m == 8;    Km=r.*sq_DE2.*(2.*E+D)./D3/sqrt(2);                                end

if m == 41;   Km=t-sq_DE2/sqrt(2);                                               end

if m == 51;   Km= -(h.*sq_DE2-t.*sq_DE1)./D/sqrt(2);                          end

if m == 61;   Km=1-(h.*sq_DE1+t.*sq_DE2)./D/sqrt(2);                          end

end

