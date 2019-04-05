clear all; clf reset;
global NumLegendreTerms
n=2; eps=1e-4;
h1=log(0.2);
h2=log(10);
h=h1:(h2-h1)/50:h2;
h=exp(h);
[n1,n2]=size(h); 
fid=fopen('umax.dat','w');
fprintf('\n');
for e=1:n2
 fprintf('Working on case %d (out of %d)',e,n2);
 [fi,psi,t,Wt]=fredholm(h(e),n,eps);
 line(t,fi), hold on
 line(t,psi), hold on
 r=0.001:0.1:2;
 [Uz,Ur]=intgr(r,fi,psi,h(e),Wt,t);
 clf
 line(r,-Uz), hold on
 line(r,Ur), hold on
 mu(e)=max(abs(Uz));
 mh(e)=max(abs(Ur));
 fprintf(fid,'%d %d %d \n',h(e),mu(e),mh(e));
 for ii=1:33, fprintf('\b'); end
end
fprintf('\n');
status=fclose(fid);
