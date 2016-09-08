function angleInRadians = deg2rad(angleInDegrees)
% DEG2RAD Convert angles from degrees to radians
%
%   angleInRadians = DEG2RAD(angleInDegrees) converts angle units from
%   degrees to radians.
%
%   See also: fromDegrees, fromRadians, toDegrees, toRadians, rad2deg.

% Copyright 2007 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2007/02/11 05:47:38 $

angleInRadians = (pi/180) * angleInDegrees;
