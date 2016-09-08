function [ellipsoid,msg] = geoidtst(ellipsoid)
%GEOIDTST  Tests for a valid ellipsoid vector
%
%   ELLIPSOID = GEOIDTST(ELLIPSOID) ensures ELLIPSOID has the form
%   [semimajor-axis-length eccentricity] where is a real non-negative
%   scalar (zero is allowed to serve as a special code) and eccentricity
%   is real scalar in the half-open interval [0 1).  If ELLIPSOID itself
%   is a scalar, then a zero eccentricity is appended.
%
%   See also private/checkellipsoid.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.10.4.2 $  $Date: 2007/10/10 20:47:42 $

% Deprecated syntax
% -----------------
% [ELLIPSOID, MSG] = GEOIDTST(ELLIPSOID) returns a string MSG indicating any
% error encountered.

% Implementation Note
% -------------------
% This is an older function that will be deprecated soon.  The
% implementation below uses the newer ignoreComplex and checkellipsoid
% functions, yet preserves the message strings used in the original
% version.

% Initialize output message string in case it's needed.
msg = '';

% Handle complex and non-numeric input
try
    ellipsoid = ignoreComplex(ellipsoid, mfilename, 'ellipsoid');
catch exception
    eid = exception.identifier;
    if strcmp(eid, 'map:geoidtst:nonNumericInput')
        exception = MException(eid, ...
            'Geoid vector must have 1 or 2 elements');
    end
    if nargout > 1
        msg = exception.message; return
    else
        exception.throw
    end
end

% Quietly convert negative semimajor axis length to positive.
% (Stay out of trouble, but avoid introducing a new error condition.)
if ellipsoid(1) < 0
    ellipsoid(1) = abs(ellipsoid(1));
end

% Finish up with CHECKELLIPSOID
try
    ellipsoid = checkellipsoid(ellipsoid, 'geoidtst', 'ellipsoid', 1);
catch exception
    % Preserve the original message text in necessary
    eid = exception.identifier;
    if strcmp(eid, 'map:geoidtst:ellipsoidNot1By2')
        exception = MException(eid, ...
            'Geoid vector must have 1 or 2 elements');
    elseif strcmp(eid, 'map:geoidtst:invalidEccentricity') || ...
            strcmp(eid, 'map:geoidtst:invalidEllipsoid')
        exception = MException(eid, ...
            'Geoid eccentricity must be in [0,1]');
    elseif strcmp(eid, 'MATLAB:geoidtst:expectedNonnegative') && ...
            ellipsoid(2) < 0
        exception = MException(eid, ...
            'Geoid eccentricity must be in [0,1]');
    end
    if nargout > 1
        msg = exception.message; return
    else
        exception.throw
    end
end
