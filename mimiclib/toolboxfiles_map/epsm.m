function epsilon = epsm(angleUnits)
%EPSM  Accuracy in angle units for certain map computations
%
%   e = EPSM returns the accuracy, in degrees, of certain computations
%   performed in the Mapping Toolbox.
%
%   e = EPSM(angleUnits) returns the accuracy in the specified angle
%   units, which can be 'degrees' or 'radians'.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.3 $  $Date: 2007/03/27 19:12:00 $

epsilon = 1.0E-6;
if nargin > 0
    epsilon = fromDegrees(angleUnits, epsilon);
end
