function [fi,psi,t,Wt]=fredholm(h,m,er)
% function [fi,psi,t,Wt]=fredholm(h,m,er)
% fi,psi: basis functions
% t: interval of integration
% m: size(t)
%er=1e-7;
lamda=2/pi;
RtWt;
NumLegendreTerms=length(Root);
t=zeros(1,m*NumLegendreTerms);
Wt=t;
for k=1:m
    for i=1:NumLegendreTerms
        d1=1/m;
        t1=d1*(k-1);
        r1=d1*k;
        j=NumLegendreTerms*(k-1)+i;
        t(j)=Root(i)*(r1-t1)*0.5+(r1+t1)*0.5;
        Wt(j)=0.5/m*Weight(i);
    end
end
fi1=-lamda*t;
psi1=zeros(size(t));
fi  =zeros(size(t));
psi =zeros(size(t));
res=1e9;

while res > er
    % To speed up, use matrix operation instead of loop, Yunjun, 2015-11-10
    Wt2  = repmat(Wt,  m*NumLegendreTerms,1);
    fi2  = repmat(fi1, m*NumLegendreTerms,1);
    psi2 = repmat(psi1,m*NumLegendreTerms,1);
    rr   = repmat(t,   m*NumLegendreTerms,1);
    tt   = repmat(t',size(t));
    fpk = fpkernel(h,tt,rr,1:4);
    fi  = sum(Wt2.*(fi2 .*fpk(:,:,1) + psi2.*fpk(:,:,3)),2)' - t;
    psi = sum(Wt2.*(psi2.*fpk(:,:,2) + fi2 .*fpk(:,:,4)),2)';
%     fi  = sum(Wt2.*(fi2.*fpkernel(h,tt,rr,1)+...,
%               psi2.*fpkernel(h,tt,rr,3)),2)' - t;
%     psi = sum(Wt2.*(psi2.*fpkernel(h,tt,rr,2)+...,
%               fi2.*fpkernel(h,tt,rr,4)),2)';
%     for i=1:m*NumLegendreTerms
%         fi(i)  = -t(i)+sum(Wt.*(fi1.*fpkernel(h,t(i),t,1)+...,
%                  psi1.*fpkernel(h,t(i),t,3)));
%         psi(i) = sum(Wt.*(psi1.*fpkernel(h,t(i),t,2)+...,
%                  fi1.*fpkernel(h,t(i),t,4)));
%     end

    fi=fi*lamda;  psi=psi*lamda;
    % find maximum relative change
    [fim,im]=max(abs(fi1-fi));
    fim=fim/abs(fi(im));
    [psim,im]=max(abs(psi1-psi));
    psim=psim/abs(psi(im));
    res=max([fim psim]);
    fi1=fi;
    psi1=psi;
end %while

%function [t,Wt] =SimpRtWt(m);  % nodes and weights for Simpson integration
%        t=0:1/m:1;
%        Wt=2/3/m*ones(size(t));
%        I=1:m+1;
%        ev=find(mod(I,2)==0);
%        Wt(ev)=4/3/m;
%        Wt(1)=1/3/m;
%        Wt(m+1)=1/3/m;

