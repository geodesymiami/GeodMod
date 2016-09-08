function output = abstractAngleConv(...
    homeUnits, alternateUnits, basicConvFcn, toOrFromUnits, varargin)
% Abstraction of operations needed by the following degrees-radians
% angle unit converters:  fromDegrees, fromRadians, toDegrees, and
% toRadians.  The first three arguments should be hard-coded within a given
% converter, which should then pass through its remaining arguments.  The
% interpretation of the first three arguments is as follows:
%
%    homeUnits       Angle units given in the function name
%                    (e.g., 'radians' for toRadians or fromRadians)
%
%    alternateUnits  Angle units not given in function name
%                    (e.g., 'degrees' for toRadians or fromRadians)
%
%    basicConvFcn    Function handle to deg2rad or rad2deg, depending
%                    on both home units and direction of conversion:
%
%                       @deg2rad for fromDegrees and toRadians
%                       @rad2deg for fromRadians and toDegrees
%
% Note that homeUnits and alternateUnits are defined relative to the
% calling function name, so they do not have a fixed correspondence to a
% more literal notion of 'fromUnits' and 'toUnits'.  As a result, this
% function may seem somewhat more abstract than necessary, but the
% reward is efficiency -- the calling functions consist of a single line
% with no string comparisons or logic of their own.

% Copyright 2007 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2007/03/27 19:14:34 $

if strncmpi(toOrFromUnits, homeUnits, numel(toOrFromUnits))
    output = varargin;
elseif strncmpi(toOrFromUnits, alternateUnits, numel(toOrFromUnits))
    output = cellfun(basicConvFcn, varargin, 'UniformOutput', false);
else
    badAngleUnits(toOrFromUnits)
end
