function [rng, az] = distance(varargin)
%DISTANCE  Distance between points on sphere or ellipsoid
%
%   [RNG, AZ] = DISTANCE(LAT1,LON1,LAT2,LON2) computes the great circle
%   distance(s) and azimuth(s) between pairs of points on the surface of
%   a sphere.  The input latitudes and longitudes, LAT1, LON1, LAT2,
%   LON2, can be scalars or arrays of equal size and must be expressed
%   in degrees.  The distance RNG is expressed in degrees of arc length
%   and will have the same size as the input arrays.  Azimuth AZ is
%   clockwise from north from the first point to the second point, in
%   units of degrees.  When given a combination of scalar and array
%   inputs, the scalar inputs are automatically expanded to match the
%   size of the arrays.
%
%   [RNG, AZ] = DISTANCE(LAT1,LON1,LAT2,LON2,ELLIPSOID) computes the
%   geodesic distance and azimuth assuming that the points lie on the
%   reference ellipsoid defined by the input ELLIPSOID.  The ellipsoid
%   vector is of the form [semimajor axis, eccentricity].  The output
%   RNG is expressed in the same distance units as the semimajor axis of
%   the ellipsoid vector.
%
%   [RNG, AZ] = DISTANCE(LAT1,LON1,LAT2,LON2,UNITS) uses the input
%   string UNITS to define the angle unit of the outputs RNG and AZ and
%   the input latitude-longitude coordinates.  UNITS may equal 'degrees'
%   (the default value) or 'radians'.
%
%   [RNG, AZ] = DISTANCE(LAT1,LON1,LAT2,LON2,ELLIPSOID,UNITS) uses the
%   UNITS string to specify the units of the latitude-longitude
%   coordinates, but the output range has the same units as the
%   semimajor axis of the ellipsoid vector.
%
%   [RNG, AZ] = DISTANCE(TRACK,...) uses the input string TRACK to
%   specify either a great circle/geodesic or a rhumb line distance
%   calculation.  If TRACK equals 'gc' (the default value), then great
%   circle distances are computed on a sphere and geodesic distances are
%   computed on an ellipsoid. If TRACK equals 'rh', then rhumb line
%   distances are computed on either a sphere or ellipsoid.
%
%   [RNG, AZ] =  = DISTANCE(PT1,PT2) accepts N-by-2 coordinate arrays
%   PT1 and PT2 such that PT1 = [LAT1 LON1] and PT2 = [LAT2 LON2] where
%   LAT1, LON1, LAT2, and LON2 are column vectors.  It is equivalent to
%   RNG = DISTANCE(PT1(:,1),PT1(:,2),PT2(:,1),PT2(:,2)).
%
%   [RNG, AZ] = DISTANCE(PT1,PT2,ELLIPSOID)
%   [RNG, AZ] = DISTANCE(PT1,PT2,UNITS),
%   [RNG, AZ] = DISTANCE(PT1,PT2,ELLIPSOID,UNITS) and
%   [RNG, AZ] = DISTANCE(TRACK,PT1,...) are all valid calling forms.
%
%   Remark on Computing Azimuths
%   ----------------------------
%   Note that when both distance and azimuth are required for the same
%   point pair(s), it's more efficient to compute both with a single
%   call to DISTANCE.  That is, use:
%
%       [RNG, AZ] = DISTANCE(...);
%
%   rather than its slower equivalent:
%
%       RNG = DISTANCE(...);
%       AZ = AZIMUTH(...);
%
%   Remark on Output Units
%   ----------------------
%   To express the output RNG as an arc length expressed in either degrees
%   or radians, omit the ELLIPSOID argument.  This is possible only on a
%   sphere.  If ELLIPSOID is supplied, RNG is a distance expressed in the
%   same units as the semimajor axis of the ellipsoid.  Specify ELLIPSOID
%   as [R 0] to compute RNG as a distance on a sphere of radius R, with
%   RNG having the same units as R.
%
%   Remark on Eccentricity
%   ----------------------
%   Geodesic distances on an ellipsoid are valid only for small
%   eccentricities typical of the Earth (e.g., 0.08 or less).
%
%   Remark on Long Geodesics
%   ------------------------
%   Distance calculations for geodesics degrade slowly with increasing
%   distance and may break down for points that are nearly antipodal,
%   and/or when both points are very close to the Equator.  In addition,
%   for calculations on an ellipsoid, there is a small but finite input
%   space, consisting of pairs in which both the points are nearly
%   antipodal AND both points fall close to (but not precisely on) the
%   Equator. In this case, a warning is issued and both RNG and AZ are
%   set to NaN for the "problem pairs."
%   
%   See also AZIMUTH, RECKON.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.10.4.8 $  $Date: 2007/12/10 21:38:26 $

[useGeodesic, lat1, lon1, lat2, lon2, ellipsoid, ...
    units, insize, useAngularDistance] = parseDistAzInputs(varargin{:});

if useGeodesic
    % Start with spherical approximation even if using an ellipsoid
    rng = greatcircledist(lat1, lon1, lat2, lon2, ellipsoid(1));
    if ellipsoid(2) == 0
        if nargout > 1
            % Azimuth on a sphere
            az = greatcircleaz(lat1, lon1, lat2, lon2);
        end
    else
        if nargout > 1
            % Distance and azimuth on an ellipsoid:
            %   Always use GEODESICINV.
            az = zeros(size(rng));
            q = (rng > 0);
            if any(q)
                [rng(q),az(q)] = geodesicinv(lat1(q), lon1(q), lat2(q), lon2(q), ellipsoid);
            end
        else
            % Distance only on an ellipsoid:
            %   Use short-geodesic approximation if possible.
            shortGeodesicLimit = 1/60;
            q = (rng/ellipsoid(1) < shortGeodesicLimit);
            if any(q)
                rng(q) = shortgeodesicdist(lat1(q), lon1(q), lat2(q), lon2(q), ellipsoid);
            end
            if any(~q)
                rng(~q) = geodesicinv(lat1(~q), lon1(~q), lat2(~q), lon2(~q), ellipsoid);
            end
        end
    end
else
    if ellipsoid(2) == 0
        [rng,az] = rhumblineinv(lat1, lon1, lat2, lon2, ellipsoid(1));
    else
        [rng,az] = rhumblineinv(lat1, lon1, lat2, lon2, ellipsoid);
    end
end

rng = reshape(rng,insize);
if useAngularDistance
    rng = fromRadians(units, rng);
end

if nargout > 1
    % Ensure azimuth in the range [0 2*pi).
    az = zero22pi(az,'radians','exact');

    % Set output units and match shape of input arrays.
    az = reshape(fromRadians(units,az),insize);
end

%----------------------------------------------------------------------

function rng = greatcircledist(lat1, lon1, lat2, lon2, r)

% Calculate great circle distance between points on a sphere using the
% Haversine Formula.  LAT1, LON1, LAT2, and LON2 are in radians.  RNG is a
% length and has the same units as the radius of the sphere, R.  (If R is
% 1, then RNG is effectively arc length in radians.)

a = sin((lat2-lat1)/2).^2 + cos(lat1) .* cos(lat2) .* sin((lon2-lon1)/2).^2;
rng = r * 2 * atan2(sqrt(a),sqrt(1 - a));

%--------------------------------------------------------------------------

function az = greatcircleaz(lat1,lon1,lat2,lon2)

% Inputs LAT1, LON1, LAT2, LON2 are in units of radians.

az = atan2(cos(lat2) .* sin(lon2-lon1),...
           cos(lat1) .* sin(lat2) - sin(lat1) .* cos(lat2) .* cos(lon2-lon1));

% Azimuths are undefined at the poles, so we choose a convention: zero at
% the north pole and pi at the south pole.
az(lat1 <= -pi/2) = 0;
az(lat2 >=  pi/2) = 0;
az(lat2 <= -pi/2) = pi;
az(lat1 >=  pi/2) = pi;

%----------------------------------------------------------------------

function rng = shortgeodesicdist(lat1, lon1, lat2, lon2, ellipsoid)

% Calculate an approximate geodesic distance between nearby points on an
% ellipsoid.  LAT1, LON1, LAT2, and LON2 are in radians.  RNG is a length
% and has the same units as the semi-major axis of the ellipsoid.

a = ellipsoid(1);              %  Semimajor axis
b = minaxis(ellipsoid);        %  Semiminor axis

% Compute the Cartesian coordinates of the points on the ellipsoid
% using their parametric latitudes.
par1 = convertlat(ellipsoid, lat1, 'geodetic', 'parametric', 'nocheck');
par2 = convertlat(ellipsoid, lat2, 'geodetic', 'parametric', 'nocheck');

[x1,y1,z1] = sph2cart(lon1,par1,a);
z1 = (b/a) * z1;

[x2,y2,z2] = sph2cart(lon2,par2,a);
z2 = (b/a) * z2;

% Compute the chord length.  Can't use norm function because
% x1, x2, etc may already be vectors or matrices.

k = sqrt( (x1-x2).^2 +  (y1-y2).^2 +  (z1-z2).^2 );

% Compute a correction factor, and then the range. The correction factor
% breaks down as the distance and/or the eccentricity increase.

r = rsphere('euler',lat1,lon1,lat2,lon2,ellipsoid,'radians');
delta = k.^3 ./ (24*r.^2) + 3*k.^5 ./ (640*r.^4);
rng = k + delta;

%--------------------------------------------------------------------------
% 
% function az = shortgeodesicaz(lat1, lon1, lat2, lon2, ellipsoid)
% 
% % Calculate great circle azimuth between points on a sphere or ellipsoid.
% %
% % Inputs LAT1, LON1, LAT2, LON2 are column vectors in units of radians.
% %
% % References
% % ----------
% % For the ellipsoid:  D. H. Maling, Coordinate Systems and
% % Map Projections, 2nd Edition Pergamon Press, 1992, pp. 74-76.
% % This forumula can be shown to be equivalent for a sphere to
% % J. P. Snyder,  "Map Projections - A Working Manual,"  US Geological
% % Survey Professional Paper 1395, US Government Printing Office,
% % Washington, DC, 1987,  pp. 29-32.
% 
% % Parametric latitudes
% cospar1 = cos(...
%     convertlat(ellipsoid, lat1, 'geodetic', 'parametric', 'nocheck'));
% cospar2 = cos(...
%     convertlat(ellipsoid, lat2, 'geodetic', 'parametric', 'nocheck'));
% 
% r = (minaxis(ellipsoid) / ellipsoid(1))^2;   % (b/a)^2
% 
% f1 = cos(lat2) .* sin(lon2-lon1);
% f2 = r * cos(lat1) .* sin(lat2);
% f3 = sin(lat1) .* cos(lat2) .* cos(lon2-lon1);
% f4 = (1 - r) * sin(lat1) .* cos(lat2) .* cospar1 ./ cospar2;
% 
% az = atan2(f1, f2 - f3 + f4);
% 
% % Azimuths are undefined at the poles, so we choose a convention: zero at
% % the north pole and pi at the south pole.
% az(lat1 <= -pi/2) = 0;
% az(lat2 >=  pi/2) = 0;
% az(lat2 <= -pi/2) = pi;
% az(lat1 >=  pi/2) = pi;
