function [latout,lonout] = reckon(varargin)
%RECKON  Point at specified azimuth, range on sphere or ellipsoid
%
%   [LATOUT, LONOUT] = RECKON(LAT, LON, RNG, AZ), for scalar inputs,
%   calculates a position (LATOUT, LONOUT) at a given range RNG and azimuth
%   AZ along a great circle from a starting point defined by LAT and LON.
%   LAT and LON are in degrees.  The range is in degrees of arc length on a
%   sphere.  The input azimuth is in degrees, measured clockwise from due
%   north.  RECKON calculates multiple positions when given four arrays of
%   matching size.  When given a combination of scalar and array inputs,
%   the scalar inputs are automatically expanded to match the size of the
%   arrays.
%
%   [LATOUT, LONOUT] = RECKON(LAT, LON, RNG, AZ, UNITS), where UNITS is any
%   valid angle units string, specifies the angular units of the inputs and
%   outputs, including RNG.  The default value is 'degrees'.
%
%   [LATOUT, LONOUT] = RECKON(LAT, LON, RNG, AZ, ELLIPSOID) calculates
%   positions along a geodesic on an ellipsoid, as specified by the
%   two-element vector ELLIPSOID.  The range, RNG, is in linear distance
%   units matching the units of the semimajor axis of the ellipsoid (the
%   first element of ELLIPSOID).
%
%   [LATOUT, LONOUT] = RECKON(LAT, LON, RNG, AZ, ELLIPSOID, UNITS)
%   calculates positions on the specified ellipsoid with LAT, LON, AZ,
%   LATOUT, and LONOUT in the specified angle units. 
%
%   [LATOUT,LONOUT] = RECKON(TRACK,...) calculates positions on great
%   circles (or geodesics) if TRACK is 'gc' and along rhumb lines if TRACK
%   is 'rh'. The default value is 'gc'.
%
%   See also AZIMUTH, DISTANCE.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.13.4.9 $  $Date: 2007/11/09 20:25:12 $

[useGeodesic,lat,lon,rng,az,ellipsoid,units,insize] = parseInputs(varargin{:});
if useGeodesic
    if ellipsoid(2) ~= 0
        [newlat,newlon] = geodesicfwd(lat,lon,az,rng,ellipsoid);
    else
        [newlat,newlon] = greatcirclefwd(lat,lon,az,rng,ellipsoid(1));
    end
else
    if ellipsoid(2) ~= 0
        [newlat,newlon] = rhumblinefwd(lat,lon,az,rng,ellipsoid);
    else
        [newlat,newlon] = rhumblinefwd(lat,lon,az,rng,ellipsoid(1));
    end
end

newlon = npi2pi(newlon,'radians');
[newlat, newlon] = fromRadians(units, newlat, newlon);
newlat = reshape(newlat, insize);
newlon = reshape(newlon, insize);

%  Set the output arguments
if nargout <= 1
    % Undocumented command-line output for single points.
    latout = [newlat newlon];
elseif nargout == 2
    latout = newlat;
    lonout = newlon;
end

%--------------------------------------------------------------------------
% 
% function lambda = wrapToPi(lambda)
% % Wrap angle in radians to [-pi pi]
% 
% gePi = (lambda >= pi);
% lambda = mod(lambda + pi, 2*pi) - pi;
% lambda(gePi & (lambda == -pi)) = pi;

%--------------------------------------------------------------------------

function [useGeodesic,lat,lon,rng,az,ellipsoid,units,insize] ...
                                              = parseInputs(varargin)

% Handle optional first input argument.
if (nargin >= 1) && ischar(varargin{1})
    trackstr = varargin{1};
    varargin(1) = [];
else
    trackstr = 'gc';
end

% Check TRACKSTR.
switch(lower(trackstr))
    case {'g','gc'}
        useGeodesic = true;
    case {'r','rh'}
        useGeodesic = false;
    otherwise
        error(['map:' mfilename ':mapError'], 'Unrecognized track string')
end

% Check argument counting after stripping off first argument.
n = numel(varargin);
error(nargchk(4,6,n,'struct'))

% Assign the fixed arguments.
lat = varargin{1};
lon = varargin{2};
rng = varargin{3};
az  = varargin{4};
                         
% Parse the optional arguments: ELLIPSOID and UNITS.
ellipsoid = [];
units     = 'degrees';

if (n == 5)
    if ischar(varargin{5})
        units = varargin{5};
    else
        ellipsoid = varargin{5};
    end
end

if (n == 6)
	ellipsoid = varargin{5};
    units     = varargin{6};
end

% If ELLIPSOID was omitted, use a unit sphere.
if isempty(ellipsoid)
    ellipsoid = [1 0];
    % Make sure RNG is a distance on the unit sphere.
    rng = toRadians(units,rng);
end

% Check for matching lat-lon-az-rng sizes and expand scalar inputs if
% necessary.
[lat, lon, az, rng, insize] = expandScalarInputs(lat, lon, az, rng);
if isempty(insize)
    error(['map:' mfilename ':inconsistentLatLonAzRngSizes'], ...
        'Lat, long, azimuth and range inputs must have same dimension')
end

% Make sure angles are in radians and convert the input
% arrays to column vectors.
[lat, lon, az] = toRadians(units, lat(:), lon(:), az(:));
rng = rng(:);

% Check the ellipsoid.
if ellipsoid(1) == 0
    % Ensure a nonzero semimajor axis (note: should be an error)
    % warning('Semimajor axis of ELLIPSOID must be positive, reset to 1.');
    ellipsoid(1) = 1;
end
ellipsoid = geoidtst(ellipsoid);
