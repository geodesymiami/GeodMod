function lon = zero22pi(lon, angleunits, method) %#ok<INUSD>
%ZERO22PI Wrap longitudes to [0 360] degree interval
%
%   ZERO22PI has been replaced by wrapTo360 (for use when working in
%   degrees) and wrapTo2Pi (for use when working in radians).
%
%   NEWLON = ZERO22PI(LON) wraps angles in degrees into the 0 to 360
%   degree range.
%
%   NEWLON = ZERO22PI(LON, ANGLEUNITS) works in the units defined by
%   the string ANGLEUNITS, which can be either 'degrees' or 'radians'.
%   ANGLEUNITS may be abbreviated and is case-insensitive.
%
%   Note
%   ----
%   ZERO22PI also accepts an optional third argument, but it is unused
%   and does not affect the resulting computation in any way.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.4 $  $Date: 2007/05/10 13:47:40 $

% Note: We determined during work on R2007b that use of the third input
% argument, METHOD, is not implemented. Instead, ZERO22PI simply passes
% the hard-code string 'exact' to NPI2PI.  However, because this
% function has been replaced, we will not expand its behavior now.
% Instead, we've removed the three-input syntax from the help.  METHOD
% remains in the function declaration, however, to ensure backward
% compatibility -- so that call of the form zero22pi(..., 'exact') would
% contiune to work without errors.

error(nargchk(1, 3, nargin, 'struct'))
if nargin == 1
    angleunits = 'degrees';
end

%  Convert inputs to radians, then wrap to the -pi to pi range
lon = toRadians(angleunits, lon);
lon = npi2pi(lon,'radians','exact');

%  Allow points near zero to remain there and shift the points in the
%  -pi to 0 range to the pi to 2pi range
epsilon = -epsm('radians');
nearZero = (lon < epsilon);
lon(nearZero) = lon(nearZero) + 2*pi;

%  Reset near-zero points
lon(lon < 0) = 0;

%  Restore original units
lon = fromRadians(angleunits, lon);
