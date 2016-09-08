function r = rsphere(varargin)
%RSPHERE  Radii of auxiliary spheres
%
%  r = RSPHERE('axes',geoid,'method') computes either the biaxial
%  or triaxial average of the semimajor and semiminor axes of
%  the ellipsoid.  If the method string is 'biaxial', then the
%  biaxial average, (a+b)/2, is returned.  If the method string is
%  'triaxial', then the average (2a+b)/3 is returned.  If the method
%  string is omitted, then the biaxial average is computed.
%
%  r = RSPHERE('eqavol',geoid) computes the radius of the sphere with
%  equal volume as the ellipsoid specified by geoid.
%
%  r = RSPHERE('authalic',geoid) computes the radius of the sphere with
%  equal surface area as the ellipsoid specified by geoid.
%
%  r = RSPHERE('rectifying',geoid) computes the radius of the sphere
%  with equal meridional distances as the ellipsoid specified by geoid.
%
%  r = RSPHERE('curve',geoid) and r = RSPHERE('curve',geoid,lat)
%  computes the arithmetic average of the transverse and meridional
%  radii of curvature at the specified latitude point.  If the input
%  lat is omitted, then lat = 45 degrees is assumed.
%
%  r = RSPHERE('curve',geoid,'method') and r = RSPHERE('curve',...
%  geoid,lat,'method) uses the averaging method specified by the
%  input string 'method' at a the specified latitude.  If lat is omitted,
%  then lat = 45 degrees is assumed.  If method = 'mean', then
%  arithmetic averaging is used.  If method = 'norm', then geometric
%  averaging is used.  r = RSPHERE('curve',geoid,lat,'method','units')
%  uses the units specified by 'units' for the input latitude.  If
%  omitted, 'degrees' are assumed.
%
%  r = RSPHERE('euler',lat1,lon1,lat2,lon2,geoid) computes the Euler
%  radii of curvature at the midpoint of the great circle arc defined
%  by the endpoints (lat1,lon1) and (lat2,lon2).
%
%  r = RSPHERE('euler',....,'units') uses the units specified by 'units'
%  for the input latitudes and longitudes.  If omitted, 'degrees' are
%  assumed.
%
%  See also RCURVE.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.11.4.4 $  $Date: 2007/11/09 20:25:20 $
% Written by:  E. Byrns, E. Brown

%  Reference: D. H. Maling, Coordinate Systems and Map
%             Projections, 2nd Edition Pergamon Press, 1992, pp.. 77-79.

error(nargchk(2, inf, nargin, 'struct'))

switch varargin{1}
case 'biaxial'           %  2 axis ellipsoid calculation
    geoid = geoidtst(varargin{2});

	method = [];    %  Initialize potential default inputs
    if nargin >= 3
        method = lower(varargin{3});
    end

    if isempty(method)            %  Test the method string
	     method = 'mean';
    else
	     strmat =['mean';'norm'];
         indx = strmatch(method,strmat);
         if length(indx) == 1;
             method = strmat(indx,:);
         end
    end

	b = minaxis(geoid);   %  Compute semiminor axis

    switch method
	    case 'mean'
            r = (geoid(1) + b) / 2;
	    case 'norm'
            r = sqrt(geoid(1) .* b);
        otherwise
            error(['map:' mfilename ':mapError'], 'Unrecognized method string');
    end

case 'triaxial'           %  3 axis ellipsoid calculation
    geoid = geoidtst(varargin{2});

	method = [];    %  Initialize potential default inputs
    if nargin >= 3
        method = lower(varargin{3});
    end

    if isempty(method)            %  Test the method string
	     method = 'mean';
    else
	     strmat =['mean';'norm'];
         indx = strmatch(method,strmat);
         if length(indx) == 1;
             method = strmat(indx,:);
         end
    end

	b = minaxis(geoid);   %  Compute semiminor axis

    switch method
	    case 'mean'
            r = (2*geoid(1) + b) / 3;
	    case 'norm'
            r = (geoid(1).^2 .* b) .^ (1/3);
        otherwise
            error(['map:' mfilename ':mapError'], 'Unrecognized method string');
    end

case 'curve'       %  Averaged radii of curvatures

    geoid = geoidtst(varargin{2});

    %  Initialize potential defaults
     
	method = [];
    lat = [];
    units = [];

   if nargin > 2
       switch nargin
           case 3
               if ischar(varargin{3})
                   method = lower(varargin{3});
               else
                   lat = varargin{3};
               end
           case 4
               lat = varargin{3};
               method = lower(varargin{4});
           otherwise
               lat = varargin{3};
               method = lower(varargin{4});
               units = varargin{5};
       end
    end

%  Empty argument tests

    if isempty(units)
        units  = 'degrees';
    end
    
    if isempty(lat)
        lat = fromDegrees(units, 45);
    end

    if isempty(method)
	     method = 'mean';
    else
	     strmat =['mean';'norm'];
         indx = strmatch(method,strmat);
         if length(indx) == 1;
             method = strmat(indx,:);
         end
    end

%  Compute the meridional and transverse radii of curvature

    rho = rcurve('meridian',geoid,lat,units);
	nu  = rcurve('transverse',geoid,lat,units);

%  Compute the radius of the sphere

    switch method
	    case 'mean',    r = (rho + nu) / 2;
	    case 'norm',    r = sqrt(rho .* nu);
		otherwise,      error(['map:' mfilename ':mapError'], 'Unrecognized method string');
   end


case 'eqavol'       %  Equal Volume Sphere
    geoid = geoidtst(varargin{2});
	flat = ecc2flat(geoid);
    r = geoid(1) * (1 - flat/3 - flat.^2 /9);

case 'authalic'       %  Equal Surface Area Sphere
    geoid = geoidtst(varargin{2});
    e = geoid(2);
    if e > 0
		fact1 = geoid(1)^2 / (2);
        fact2 = (1 - e^2) / (2*e);
        fact3 = log((1+e) / (1-e));
        r = sqrt(fact1 * (1 + fact2 * fact3));
	else
		r = geoid(1);
    end

case 'rectifying'       %  Equal Meridian Distance Sphere
    geoid = geoidtst(varargin{2});
    n = ecc2n(geoid);
    fact1 = (1 - n) * (1 - n^2);
    r = geoid(1) * fact1 * (1 + 9*n^2 /4 + 225*n^4 /64);

case 'euler'        %  Euler radius of curvature

    if nargin < 6
        error(['map:' mfilename ':mapError'], ...
            'Incorrect number of arguments')
    end

    geoid = geoidtst(varargin{6});

    if nargin >= 7
        units = varargin{7};
    else
        units = 'degrees';
    end

%  Convert inputs to radians

    [lat1, lon1, lat2, lon2] = toRadians(units, varargin{2:5});

%  Compute the mid-latitude point

    latmid = lat1 + (lat2 - lat1)/2;

%  Compute the azimuth

    az = azimuth('gc',lat1,lon1,lat2,lon2,geoid,'radians');

%  Compute the meridional and transverse radii of curvature

    rho = rcurve('meridian',geoid,latmid,'radians');
	nu  = rcurve('transverse',geoid,latmid,'radians');

%  Compute the radius of the arc from point 1 to point 2
%  Ref:  Maling, p. 76.

   den = rho .* sin(az).^2 + nu .* cos(az).^2;
   r = rho .* nu ./ den;

otherwise
   error(['map:' mfilename ':mapError'], 'Unrecognized calculation string')

end
