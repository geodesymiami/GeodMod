function [G1, G2, G3] = MakeKernel(patch_geom, slip, xy, nu)

%% [G1, G2, G3] = MakeKernel(patch_geom, slip, xy, nu)

%Subroutine to make design matrix for slip

%Input:
%       patch_geom      - Patch Geometry
%	slip		- [strike-slip, dip-slip, opening]
%       xy              - xy coordinates of statinos
%       nu              - Poisson's ratio

        [nf,nel] = size(patch_geom);
        [nsta,nel2] = size(xy);
        ndata = 3*nsta; 
	
	G1 = []; G2 = []; G3 = [];          

	if slip(1)~=0
		G1 = zeros(ndata,nf);
		for i=1:nf
         		u_rel_1 = abs_disp(nu, [patch_geom(i,:),1,0,0], xy);
         		G1(:,i) = u_rel_1(:);
       		end
	end

	if slip(2)~=0
		G2 = zeros(ndata,nf);
		for i=1:nf
         		u_rel_2 = abs_disp(nu, [patch_geom(i,:),0,1,0], xy);
          		G2(:,i) = u_rel_2(:);
       	 	end
	end

	if slip(3)~=0
		G3 = zeros(ndata,nf);
		for i=1:nf
         		u_rel_3 = abs_disp(nu, [patch_geom(i,:),0,0,1], xy);
         		G3(:,i) = u_rel_3(:);
       		 end
	end


