function [radarlookvec]=GenerateLOSVec(incanglevec,azimuth)
%GenerateLOSvector      - returns LOS field based on corner co-ordinates
%
%usage:
%
%Output:  radarlookvec     [nx3] LOS vector (ENU) for n co-ordinates
%         incanglevec      [nx1] incidence angle vector for n co-ordinates
%    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radarlookvec = [-sind(incanglevec(:))*cosd(azimuth) sind(incanglevec(:))*sind(azimuth) cosd(incanglevec(:))]';


