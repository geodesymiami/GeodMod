function varargout = fromRadians(toUnits, varargin)
%fromRadians Convert angles from radians
%
%   [angle1, angle2, ...]
%       = fromRadians(toUnits, angle1InRadians, angle2InRadians, ...)
%   converts angle1InRadians, angle2InRadians, ... from radians to the
%   specified output ("to") angle units.  toUnits can be either
%   'degrees' or 'radians' and may be abbreviated.  The inputs
%   angle1InRadians, angle2InRadians, ... and their corresponding
%   outputs are numeric arrays of various sizes, with size(angleN)
%   matching size(angleNInRadians).
%
%   See also: fromDegrees, rad2deg, toDegrees, toRadians.

% Copyright 2007 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2007/02/11 05:47:40 $

varargout = abstractAngleConv( ...
    'radians', 'degrees', @rad2deg, toUnits, varargin{:});
