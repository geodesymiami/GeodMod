function lon = npi2pi(lon, angleunits, method)
%NPI2PI Wrap latitudes to [-180 180] degree interval
%
%   NPI2PI has been replaced by wrapTo180 (for use when working in
%   degrees) and wrapToPi (for use when working in radians).
%
%   NEWLON = NPI2PI(LON) wraps angles in degrees into the
%   -180 to 180 degree range.
%
%   NEWLON = NPI2PI(LON, ANGLEUNITS) works in the units defined by the
%   string ANGLEUNITS, which can be either 'degrees' or 'radians'.
%   ANGLEUNITS may be abbreviated and is case-insensitive.
%
%   NEWLON = NPI2PI(LON, ANGLEUNITS, METHOD) allows special alternative
%   computations to be used when NPI2PI is called from within certain
%   Mapping Toolbox functions.   METHOD can be one of the following
%   strings: 'exact' for exact wrapping (the default value);  'inward',
%   where angles are scaled by a factor of (1 - epsm('radians')) before
%   wrapping;  and 'outward', where angles are scaled by a factor of
%   (1 + epsm('radians')) before wrapping.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.4 $  $Date: 2007/05/10 13:47:35 $

% A remark in older versions of the reference pages suggests that the
% 'inward' and 'outward' methods were introduced purely as an attempt at
% optimization.
%
% Also, to paraphrase notes from the original authors:
%
%   If METHOD is 'inward' or 'outward', ATAN2 is used.  This may map
%   -pi, -3*pi, -5*pi, ... to pi rather than -pi because, depending on
%   the platform, atan2(-1,0) may evaluate to either -pi or pi.

error(nargchk(1, 3, nargin, 'struct'))
if nargin == 1
    angleunits = 'degrees';
    method = 'exact';
elseif nargin == 2
    method = 'exact';
end

%  Convert inputs to radians
lon = toRadians(angleunits,lon);

switch lower(method)
case 'exact'

%  The 'exact' method is not straightforward because we require that
%  odd, positive multiples of pi map to pi and odd, positive multiples
%  of -pi map to -pi.  For example:
%
%  -3pi maps to -pi;  -pi maps to -pi;  pi maps to pi;  3pi maps to pi
%
%  (Actually, this is easy to do with MOD, as in wrapToPi and wrapTo180,
%  but there a subtle numerical differences.  Thus the following code
%  remains here so as not to alter existing behavior in any way.)

    lon = pi*((abs(lon)/pi) - 2*ceil(((abs(lon)/pi)-1)/2)) .* sign(lon);

case 'inward'

%  Move values toward zero, then call atan2.

    lon = (1 - epsm('radians')) * lon;
	lon = atan2(sin(lon),cos(lon));

case 'outward'

%  Move values away from zero, then call atan2.

    lon = (1 + epsm('radians')) * lon;
	lon = atan2(sin(lon),cos(lon));

otherwise
    error('map:npi2pi:unknownMethod', ...
          'Unrecognized METHOD string: %s.', method)
end

%  Convert back to original units
lon = fromRadians(angleunits,lon);  
