% Sjonni's old code.

plot=true;
if ~plot return ; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------
% Script to plot up the estimated slip distribution as patches in 3D
%-------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load colormap
load slipcol

%plot fault geometry
[fx,fy,fz]=flakes(pm');
%figure; fill3(fx,fy,-fz,[1:size(pm,1)]); axis image

nf=size(pm,1);
sslip = s(1:nf);
if slip(1)==1
   dslip=s(nel+1:2*nel);
else
   dslip=zeros(size(sslip));
end
   
%--------------------------------------------------------
% if Dip-slip plot separetely
%if length(s)>nf+12
if slip(1)==1      
      % Plot Strike slip and dip slip distribution
      %dslip = s(nf+1:2*nf);
      figure
      subplot(211)
      fill3(fx,fy,-fz,sslip'); colormap(slipcol); colorbar; axis image
      title(['Predicted Strike-slip Distribution, roughness = ', num2str(100*roughness),' cm/km']); hold on
      %plot3(hypo(:,1),hypo(:,2),hypo(:,3),'*k')	

      subplot(212)
      fill3(fx,fy,-fz,dslip'); colormap(slipcol); colorbar; axis image
      title(['Predicted Dip-slip Distribution, roughness = ', num2str(100*roughness),' cm/km']); hold on
      %plot3(hypo(:,1),hypo(:,2),hypo(:,3),'*k')	
      
      % Plot Slip magnitude and Rake

      slipmag  = (sslip.^2 + dslip.^2).^(1/2);
      sliprake = 180 - atan(dslip./(sslip+1e-6))*180/pi;
      
      figure
      subplot(211)
      fill3(fx,fy,-fz,slipmag'); colormap(slipcol); colorbar; axis image
      title('Predicted Slip Magnitude'); hold on
      %plot3(hypo(:,1),hypo(:,2),hypo(:,3),'*k')	

      subplot(212)
      fill3(fx,fy,-fz,sliprake'); colormap(slipcol); colorbar; axis image
      title('Predicted Rake'); hold on
      %plot3(hypo(:,1),hypo(:,2),hypo(:,3),'*k')	

    
else  % IF strike slip only
      figure      
      fill3(fx,fy,-fz,s'); colormap(slipcol); colorbar; axis image;
      title(['Predicted Strike-slip Distribution, Roughness = ',num2str(roughness*100), ' cm/km']); hold on
      %plot3(hypo(:,1),hypo(:,2),hypo(:,3),'*k')	
end


nn=[nel]; ca = [0 max(abs(s))];
if slip(3)==1;
   [ffx,ffz]=flakes2D(pm',nhorz,nvert);
   figure;set(gcf,'Position',[428 68 460 560]); colormap(slipcol)
   subplot(211); patch(ffx,ffz,s'); axis image; hold on; caxis(ca);
   title('North     (strike-slip)   South'); colorbar;ax=axis;%plot(hypo(1,2),hypo(1,3),'*');ax=axis;
   ylabel('Depth down-dip (km)'); xlabel('Distance along fault')
end

