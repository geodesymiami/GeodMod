j=8;
i=63;
load S1
      [U,D,sigma,flag1]=disloc3d(PM(:,i),x,mu,0.25);
      t=[nhat'*sigma([1 2 3],:)
         nhat'*sigma([2 4 5],:)
         nhat'*sigma([3 5 6],:)];
      S1(:,i)=(nhat'*t)';
      pm1=PM;
      x1=x;
      sigma1=sigma;
t1=t;
load S1
      [U,D,sigma,flag2]=disloc3d(PM(:,i),x,mu,0.25);
      t=[nhat'*sigma([1 2 3],:)
         nhat'*sigma([2 4 5],:)
         nhat'*sigma([3 5 6],:)];
      S2(:,i)=(nhat'*t)';
      pm2=PM;
      x2=x;
      sigma2=sigma;
t2=t;
find((pm1-pm2)~=0)'
find((x1-x2)~=0)'
find((flag1-flag2)~=0)'
find((sigma1-sigma2)~=0)'
[sigma1(1), sigma2(1)]

