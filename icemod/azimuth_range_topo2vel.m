function [Udisp,mask,aspect,slope] = azimuth_range_topo2vel(dem,Uazi,Ulos,icethick,cutoffang)

%
%   vel = demazilos2vel(dem,mai,los,icethick,cutoff_hor,cutoffver)
%
% Compute velocity assuming surface parallel flow
%
%  Inputs are:
%   dem:       Topography
%   Uazi:       Azimuth motion
%   Ulos:       Line of sight motion
%   icethick:  Ice thickness (purpose of smoothing the topography)
%   cutoffang: Cutoff angles between displacement and look vectors 
%
%  Outputs are:
%   Udisp          : Ice velocity assuming parralel flow
%   mask           : Mask corresponding to the orthogonal direction from the plan azimuth/los
%   slopeDirection : Direction of greater terrain slope (0-360 clockwise from north)
%   slope          : Terrain slope (positives from the horizontal)
%
%
% N. Gourmelen - June 2009
%


%% Smooth the dem over 10 ice thickness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(dem.x_unit,'degres')
    lat_dist    = distance(dem.y_first,dem.x_first,dem.y_first+dem.y_step,dem.x_first,almanac('earth','WGS84')) ;
    lon_dist    = distance(dem.y_first,dem.x_first,dem.y_first,dem.x_first+dem.x_step,almanac('earth','WGS84')) ;
    dist        = (lat_dist+lon_dist) / 2 * 1000 ;
    corr_length = 10*icethick / dist ;
end

dem_smooth = single(lowpassfilter(dem.data,corr_length));

%% Compute DEM slope parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[LONm,LATm] = meshgrid (single(dem.x_first:dem.x_step:dem.x_first+(dem.width-1)*dem.x_step),single(dem.y_first:dem.y_step:dem.y_first+(dem.file_length-1)*dem.y_step)) ;
[aspect, slope, gradS, gradW] = gradientm(LATm,LONm,dem_smooth) ;

% Rotate N, E into azimuth, range reference frame
%rotang = 180 - Ulos.heading ;                   %% Works only for descending!
rotang = Ulos.heading ; % Opposite of normal sens because of Matlab's top-down reference frame
gradR = gradW*cosd(rotang)-gradS*sind(rotang) ;
gradA = gradS*cosd(rotang)+gradW*sind(rotang) ;

% Build mask hortogonal to 
mask = ones(size(slope)) ;
mask(find(slope >= 20 & slope <= 24 & aspect > mod(Ulos.heading -5 - 90,360) & aspect < mod(Ulos.heading +5 - 90,360))) = 0 ;


% Transform into azimuth and range slopes in degrees

slopeR = atand(gradR) ;
slopeA = atand(gradA) ;


%% Compute velocity

Udisp = Ulos ;
switch(Ulos.sat)
    case{'ErsD','EnvD'}
        %Udisp.data = sqrt((Uazi.data./cosd(slopeA)).^2+(Ulos.data./cosd(-slopeR+67)).^2) ;
        Udisp.data = sqrt((Uazi.data./cosd(slopeA)).^2+(Ulos.data./sind(23+slopeR)).^2) ;
    case{'ErsA','EnvA'}
        Udisp.data = sqrt((Uazi.data./cosd(slopeA)).^2+(Ulos.data./sind(23-slopeR)).^2) ;
end



