function r = rcurve(str,geoid,lat,units)
%RCURVE  Radii of curvature of ellipsoid
%
%  r = RCURVE(geoid,lat) computes the parallel radius of curvature
%  for an ellipsoid.  The parallel radius of curvature is the radius
%  of the small circle encompassing the globe at the specified latitude.
%  The input geoid is a standard 2 element geoid vector.  The returned
%  radius r is in the units of the first element of a standard geoid vector.
%
%  r = RCURVE(geoid,lat,'units')  defines the 'units' of the input
%  latitude data.  If omitted, 'degrees' are assumed.
%
%  r = RCURVE('parallel',...) also computes the parallel radius of
%  curvature for the ellipsoid.
%
%  r = RCURVE('meridian',...) computes the meridianal radius of
%  curvature for the ellipsoid.  The meridianal radius is the
%  radius of curvature at the specified latitude for the ellipse
%  defined a meridian.
%
%  r = RCURVE('transverse',...) computes the transverse radius of
%  curvature for the ellipsoid.  The transverse radius is the radius
%  of the curve formed by a plane intersecting the ellipsoid at the
%  latitude which is normal to the surface of the ellipsoid.
%
%  See also RSPHERE.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.13.4.4 $  $Date: 2007/11/09 20:25:11 $
% Written by:  E. Byrns, E. Brown

%  Reference: D. H. Maling, Coordinate Systems and Map Projections, 2nd Edition
%  Pergamon Press, 1992, p. 69.

if nargin < 2 || (nargin == 2 && ischar(str))
    error(['map:' mfilename ':mapError'], 'Incorrect number of arguments')

elseif (nargin == 2 && ~ischar(str)) || (nargin == 3 && ischar(str))
    if ~ischar(str)       %  Shift inputs since str omitted by user
        lat = geoid;     geoid = str;     str = [];
    end

	units = [];

elseif nargin == 3 && ~ischar(str)

    %  Shift inputs since str omitted by user
    units = lat;   lat = geoid;     geoid = str;     str = [];
end


%  Empty argument tests

if isempty(str)
     str   = 'parallel';
else
     validstrs = {'parallel','meridian','transverse'};
     indx = strmatch(lower(str),validstrs);
     if length(indx) == 1
	      str = validstrs{indx};
	 else
	      error(['map:' mfilename ':mapError'], ...
              'Unrecognized curve string')
     end
end



%  Input tests

if isempty(units);    units = 'degrees';    end
geoid = geoidtst(geoid);

if  ndims(lat) > 2
   error(['map:' mfilename ':mapError'], ...
       'Latitude limited to 2 dimensions')
end

%  Get the semimajor axis and eccentricity

semimajor = geoid(1);
eccent    = geoid(2);

%  Transform the latitude into radians and ensure real input

lat = toRadians(units,real(lat));


switch str
case 'parallel'                  %  Parallel radius of curvature

    num = 1-eccent^2;                        %  Compute the distance from
    den = 1 - (eccent * cos(lat)).^2;        %  the center of the geoid to
    rho = semimajor * sqrt(num ./ den);      %  the specified point

    r = rho .* cos(lat);

case 'meridian'                    %  Meridional radius of curvature

    num = semimajor * (1-eccent^2);
    den = 1 - (eccent * sin(lat)).^2;
	r   = num ./ sqrt(den.^3);

case 'transverse'                  %  Transverse radius of curvature

    den = 1 - (eccent * sin(lat)).^2;
	r   = semimajor ./ sqrt(den);

end
