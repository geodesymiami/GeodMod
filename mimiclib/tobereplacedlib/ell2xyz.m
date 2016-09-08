function [x,y,z]=ell2xyz(lat,lon,h,a,e2)
% ELL2XYZ  Converts ellipsoidal coordinates to cartesian.
%   Vectorized.
% Version: 18 Jan 96
% Useage:  [x,y,z]=ell2xyz(lat,lon,h,a,e2)
% Input:   lat - vector of ellipsoidal latitudes (radians)
%          lon - vector of ellipsoidal E longitudes (radians)
%          h   - vector of ellipsoidal heights (m)
%          a   - ref. ellipsoid major semi-axis (m)
%          e2  - ref. ellipsoid eccentricity squared
% Output:  x \
%          y  > vectors of cartesian coordinates in CT system (m)
%          z /

v = a./sqrt(1-e2*sin(lat).*sin(lat)) ;
x = (v+h).*cos(lat).*cos(lon)        ;
y = (v+h).*cos(lat).*sin(lon)        ;
z = (v.*(1-e2)+h).*sin(lat)          ;
