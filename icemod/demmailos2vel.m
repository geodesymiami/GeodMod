function vel=demmailos2vel(dem,mai,los,icethick)


% Compute velocity assuming surface parallel flow
%
% N. Gourmelen
%

igram=LoadData('/RAID4/ngourmelen/Iceland/Langjokull_des_ErsDifgram/geo_adf05032_940222-940225-sim_SIM_2rlks_c15_Bi.unw',struct('origin','roi_pac'));
mai=LoadData('/RAID4/ngourmelen/Iceland/Langjokull_ErsD_mai/geo_removeiono_adf_mai_adf_940222-940225-sim_SIM_4rlks.hgt',struct('origin','roi_pac'));
dem=LoadData('/RAID4/ngourmelen/Iceland/Langjokull_des_ErsDifgram/Langjokull_mosaic2_ll.dem',struct('origin','dem'));

mai.data   = mai.data/2/pi*1000  ;
igram.data = igram.data/2/pi*2.8 ;
icethick   = 0.750 ;

ang_cutoff = 5 ;

%% Smooth the dem over 10 ice thickness

if strcmp(dem.PROJECTION,'LATLON')
lat_dist = distance(dem.y_first,dem.x_first,dem.y_first+dem.y_step,dem.x_first,almanac('earth','WGS84')) ;
lon_dist = distance(dem.y_first,dem.x_first,dem.y_first,dem.x_first+dem.x_step,almanac('earth','WGS84')) ;
dist = (lat_dist+lon_dist) / 2 ;
corr_length = 10*icethick / dist ;


dem_smooth = lowpassfilter(dem.data,corr_length);

%% Compute DEM slope parameters

[LONm,LATm] = meshgrid (dem.x_first:dem.x_step:dem.x_first+(dem.width-1)*dem.x_step,dem.y_first:dem.y_step:dem.y_first+(dem.file_length-1)*dem.y_step) ;

[aspect, slope, gradN, gradE] = gradientm(LATm,LONm,dem_smooth) ;

%% Compute velocities
%%%%%%%%%%%%%%%%%%%%%
% angle = mod(atan2(a(1)*b(2)-b(1)*a(2),a(1)*b(1)+a(2)*b(2)),2*pi)*180/pi

slope = abs(slope) ;  def_head = mod(atan2(-gradE,-gradN),2*pi)*180/pi ;

% Prepare mask for LOS velocity

los_cutoff = 0.3 ;  
los_ang_hor = (sind(mod(def_head-194,360))) ;
los_ang_ver = cosd(slope) ;

los_mask_hor = find(abs(los_ang_hor) < los_cutoff) ;
los_mask_ver = find(abs(los_ang_ver) < los_cutoff) ;

los_mask = unique([los_mask_hor ; los_mask_ver]) ;

% Prepare mask for MAI velocity

mai_cutoff = 0.3 ;  
mai_ang_hor = cosd(mod(def_head-194,360)) 	  ;  
mai_ang_ver = cosd(23+slope)                      ;
mai_mask_hor = find(abs(mai_ang_hor) < mai_cutoff)  ;
mai_mask_ver = find(abs(mai_ang_ver) < mai_cutoff)  ;
mai_mask = unique([mai_mask_hor;mai_mask_ver]) ;


los_vel = igram.data./(los_ang.*cosd(23+slope)) ;  los_vel(los_mask) = 0 ;
mai_vel = mai.data./(mai_ang.*cosd(slope)) ;  mai_vel(mai_mask) = 0 ;

vel = los_vel + mai_vel ;

vel = mai.data./(mai_ang.*cosd(slope))+igram.data./(los_ang.*cosd(23+slope)) ;

