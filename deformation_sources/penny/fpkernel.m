function [K]=fpkernel(h,t,r,n)
% Kernels calculation
% Support matrix 't','r' and array 'n' input, Yunjun, 2015-11-10

p=4*h^2;
% K=[];
%[dumb,nr]=size(r);
%[dumb,nt]=size(t);


%Support array 'n' input, Yunjun, 2015-11-10
[s1,s2] = size(t);
s3 = length(n);
K=zeros(s1,s2,s3);
kg = zeros(s1,s2,3);
kg(:,:,1) = KG( t-r, p);
kg(:,:,2) = KG( t+r, p);
kg(:,:,3) = KG(-t+r, p);
for i=1:s3
    K(:,:,i) = Kn(h,t,r,n(i),p,kg);
end
if s3 == 1;  K = reshape(K,size(t)); end     % same output as Q if len(n)==1

% switch n
%     case 1	%KN
%         K=p*h*(KG(t-r,p)-KG(t+r,p));
%         
%     case 2	%KN1
%         Dlt=1e-6;
%         a=t+r;
%         b=t-r;
%         y=a.^2;
%         z=b.^2;
%         g=2*p*h.*(p^2+6*p.*(t.^2+r.^2)+5.*(a.*b).^2);
%         s=((p+z).*(p+y)).^2;
%         s=g./s;
%         trbl=-4*h./(p+t.^2).*ones(size(t));
%         rs=find(r>Dlt);
%         if t<Dlt
%             trbl=-4*h./(p+r.^2);
%         else
%             trbl(rs)=h./t(rs)./r(rs).*log((p+z(rs))./(p+y(rs)));
%         end
%         K=trbl+s+h.*(KERN(b,p)+KERN(a,p));
%         
%     case 3	%KM
%         
%         y=(t+r).^2; z=(t-r).^2;
%         a=((p+y).*(p+z)).^2;
%         c=t+r;  d=t-r;
%         b=p.*t.*((3*p^2-(c.*d).^2+2*p.*(t.^2+r.^2))./a);
%         a=p/2*(c.*KG(c,p)+d.*KG(d,p));
%         K=b-a;
%         
%     case 4	%%KM1(t,r)=KM(r,t);
%         
%         y=(t+r).^2; z=(t-r).^2;
%         a=((p+y).*(p+z)).^2;
%         c=t+r; d=-t+r;
%         b=p.*r.*((3*p^2-(c.*d).^2+2*p.*(t.^2+r.^2))./a);
%         a=p/2*(c.*KG(c,p)+d.*KG(d,p));
%         K=b-a;
%         
% end %switch n
end

%% Support array 'n' input, Yunjun, 2015-11-10
function [Km] = Kn(h,t,r,m,p,kg)

if m == 1
    Km=p*h*(kg(:,:,1)-kg(:,:,2));
end

if m == 2
    Dlt=1e-6;
    a=t+r;      b=t-r;
    y=a.^2;     z=b.^2;
    g=2*p*h.*(p^2+6*p.*(t.^2+r.^2)+5.*(a.*b).^2);
    s=((p+z).*(p+y)).^2;
    s=g./s;
    trbl=-4*h./(p+t.^2).*ones(size(t));
    rs=find(r>Dlt);
    if t<Dlt;   trbl=-4*h./(p+r.^2);
    else        trbl(rs)=h./t(rs)./r(rs).*log((p+z(rs))./(p+y(rs)));
    end
    Km=trbl+s+h.*(KERN(b,p)+KERN(a,p));
end

if m == 3
    y=(t+r).^2; z=(t-r).^2;
    a=((p+y).*(p+z)).^2;
    c=t+r;  d=t-r;
    b=p.*t.*((3*p^2-(c.*d).^2+2*p.*(t.^2+r.^2))./a);
    a=p/2*(c.*kg(:,:,2)+d.*kg(:,:,1));
    Km=b-a;
end

if m == 4
    y=(t+r).^2; z=(t-r).^2;
    a=((p+y).*(p+z)).^2;
    c=t+r; d=-t+r;
    b=p.*r.*((3*p^2-(c.*d).^2+2*p.*(t.^2+r.^2))./a);
    a=p/2*(c.*kg(:,:,2)+d.*kg(:,:,3));
    Km=b-a;
end

end

%%
function [fKG]=KG(s,p)
z=s.^2;
y=p+z;
fKG=(3*p-z)./y.^3;
end

function [fKERN]=KERN(w,p)
u=(p+w.^2).^3;
fKERN=2*(p^2/2+w.^4-p*w.^2/2)./u;
end


