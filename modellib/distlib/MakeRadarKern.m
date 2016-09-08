function [K1,K2,K3] = MakeRadarKern(pm,slip,coord,nu,nhe,nve,radarlook,Nsites)

% FA 07/2008: modified so that it works with radarsat given as 
% unit vector (1x3) (Sjonni's old version) and as LOSfield (3xN)
%
% Sjonni Jonsson, modfied by FA

if length(radarlook)==3
   LOSfield_flag=false;
else
   LOSfield_flag=true;
end

[G1, G2, G3] = MakeKernel(pm, slip, coord', nu);

K1=[];K2=[];K3=[];

if slip(1)~=0
  K1 = zeros(Nsites,nve*nhe);
  for i =1:Nsites
      if LOSfield_flag
         K1(i,:) = radarlook(i,1:3)*G1(3*(i-1)+1:3*(i-1)+3,:);
      else
         K1(i,:) = radarlook'*G1(3*(i-1)+1:3*(i-1)+3,:);
     end
  end
end

if slip(2)~=0
  K2 = zeros(Nsites,nve*nhe);
  for i =1:Nsites
      if LOSfield_flag
         K2(i,:) = radarlook(i,1:3)*G2(3*(i-1)+1:3*(i-1)+3,:);
      else
         K2(i,:) = radarlook'*G2(3*(i-1)+1:3*(i-1)+3,:);
      end
  end
end

if slip(3)~=0
  K3 = zeros(Nsites,nve*nhe);
  for i =1:Nsites
      if LOSfield_flag
         K3(i,:) = radarlook(i,1:3)*G3(3*(i-1)+1:3*(i-1)+3,:);
      else
         K3(i,:) = radarlook'*G3(3*(i-1)+1:3*(i-1)+3,:);
      end
  end
end
