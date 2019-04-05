cc
mstar(1:7)=[10 ; 2 ; 3 ; -65 ; 230 ; 6 ; 8 ]';
mu=3e9; val=[0.25 ; mu ; 30 ; 15];
mstar(8)=-3e6;      % 3MPa

load slipcol
figure('DefaultPatchLineStyle','none') ;
coord=zeros(2,1);CamPos=[10 , -10 , 0];

mstar(9)=0e6;      % 1MPa/km
h=subplot(3,1,1)
              mlinpd=[mstar(1:9)];  mlinpd([1 2 3 6 7])=mlinpd([1 2 3 6 7])*1000; mlinpd(9)=mlinpd(9)/1000;
      [u,pm]=LinearPressureDike(mlinpd,coord*1000,val(1),val(2),val(3),val(4));
              pm([1 2 3 6 7],:)=pm([1 2 3 6 7],:)/1000;         %to return in km
      op=pm(10,:);       % opening vector
      [fx,fy,fz]=flakes(pm(1:7,:));
      str=sprintf('P=%4.1f+ %4.1f/km MPa, mu=%3.1fGPa %4.1f*%4.1fkm^2, ',mstar(8)/10^6,mstar(9)/10^6,mu/10^9,mstar(1),mstar(2)); display(str)
      fill3(fx,fy,-fz,op); colormap(slipcol); colorbar; axis image; grid on
      title(str); hold on
      set(h,'CameraPosition',CamPos)

mstar(9)=-1e6;      % 1MPa/km
h=subplot(3,1,2)
              mlinpd=[mstar(1:9)];  mlinpd([1 2 3 6 7])=mlinpd([1 2 3 6 7])*1000; mlinpd(9)=mlinpd(9)/1000;
      [u,pm]=LinearPressureDike(mlinpd,coord*1000,val(1),val(2),val(3),val(4));
              pm([1 2 3 6 7],:)=pm([1 2 3 6 7],:)/1000;         %to return in km
      op=pm(10,:);       % opening vector
      [fx,fy,fz]=flakes(pm(1:7,:));
      str=sprintf('P=%4.1f+ %4.1f/km MPa, mu=%3.1fGPa %4.1f*%4.1fkm^2, ',mstar(8)/10^6,mstar(9)/10^6,mu/10^9,mstar(1),mstar(2)); display(str)
      fill3(fx,fy,-fz,op); colormap(slipcol); colorbar; axis image; grid on
      title(str); hold on
      set(h,'CameraPosition',CamPos)

mstar(9)=1e6;      % 1MPa/km
h=subplot(3,1,3)
              mlinpd=[mstar(1:9)];  mlinpd([1 2 3 6 7])=mlinpd([1 2 3 6 7])*1000; mlinpd(9)=mlinpd(9)/1000;
      [u,pm]=LinearPressureDike(mlinpd,coord*1000,val(1),val(2),val(3),val(4));
              pm([1 2 3 6 7],:)=pm([1 2 3 6 7],:)/1000;         %to return in km
      op=pm(10,:);       % opening vector
      [fx,fy,fz]=flakes(pm(1:7,:));
      str=sprintf('P=%4.1f+ %4.1f/km MPa, mu=%3.1fGPa %4.1f*%4.1fkm^2, ',mstar(8)/10^6,mstar(9)/10^6,mu/10^9,mstar(1),mstar(2)); display(str)
      fill3(fx,fy,-fz,op); colormap(slipcol); colorbar; axis image; grid on
      title(str); hold on
      set(h,'CameraPosition',CamPos)

