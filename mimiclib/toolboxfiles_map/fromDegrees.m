function varargout = fromDegrees(toUnits, varargin)
%fromDegrees Convert angles from degrees
%
%   [angle1, angle2, ...]
%       = fromDegrees(toUnits, angle1InDegrees, angle2InDegrees, ...)
%   converts angle1InDegrees, angle2InDegrees, ... from degrees to the
%   specified output ("to") angle units.  toUnits can be either
%   'degrees' or 'radians' and may be abbreviated.  The inputs
%   angle1InDegrees, angle2InDegrees, ... and their corresponding
%   outputs are numeric arrays of various sizes, with size(angleN)
%   matching size(angleNInDegrees).
%
%   See also: deg2rad, fromRadians, toDegrees, toRadians.

% Copyright 2007 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2007/02/11 05:47:39 $

varargout = abstractAngleConv( ...
    'degrees', 'radians', @deg2rad, toUnits, varargin{:});
