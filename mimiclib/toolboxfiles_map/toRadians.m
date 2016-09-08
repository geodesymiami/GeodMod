function varargout = toRadians(fromUnits, varargin)
%toRadians Convert angles to radians
%
%   [angle1InRadians, angle2InRadians, ...]
%       = toRadians(fromUnits, angle1, angle2, ...)
%   converts angle1, angle2, ... to radians from the specified input
%   ("from") angle units.  fromUnits can be either 'degrees' or 'radians'
%   and may be abbreviated.  The inputs angle1, angle2, ... and their
%   corresponding outputs are numeric arrays of various sizes, with
%   size(angleNInRadians) matching size(angleN).
%
%   See also: deg2rad, fromDegrees, fromRadians, toDegrees.

% Copyright 2007 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2007/02/11 05:47:43 $

varargout = abstractAngleConv( ...
    'radians', 'degrees', @deg2rad, fromUnits, varargin{:});
