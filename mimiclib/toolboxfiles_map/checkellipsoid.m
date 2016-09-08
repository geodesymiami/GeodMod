function ellipsoid = checkellipsoid(ellipsoid, func_name, var_name, arg_pos)
%CHECKELLIPSOID Check validity of reference ellipsoid vector
%
%     This function is intentionally undocumented and is intended for
%     use only by other Mapping Toolbox functions.  Its behavior may
%     change, or the function itself may be removed, in a future
%     release.
%
%   ELLIPSOID = CHECKELLIPSOID(ELLIPSOID, FUNC_NAME, VAR_NAME, ARG_POS)
%   ensures that ELLIPSOID is a 1-by-2 vector of the form
%
%               [semimajor-axis-length eccentricity]
%
%   with 0 <= Eccentricity < 1.  (A scalar input is interpreted as the
%   semimajor axis length and a zero eccentricity is appended in this case.)

% Copyright 2007 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2007/10/10 20:50:12 $

validateattributes( ...
    ellipsoid, {'double'}, {'real', 'finite', 'nonnegative'}, ...
    func_name, var_name, arg_pos)

if numel(ellipsoid) == 1
    % Append zero eccentricity given scalar input
    ellipsoid(1,2) = 0;
else
    % Ensure a 1-by-2 vector
    assert(isequal(size(ellipsoid),[1 2]), ...
        ['map:' func_name ':ellipsoidNot1By2'], ...
        ['Function %s expected its %s input argument, %s,\n', ...
        'to be 1-by-1 or 1-by-2.'], ...
        upper(func_name), num2ordinal(arg_pos), var_name);
end

% Ensure eccentricity < 1
assert(ellipsoid(2) < 1, ...
    ['map:' func_name ':invalidEccentricity'], ...
    ['Function %s expected the second element (the eccentricity)\n', ...
    'of %s to be in the range [0 1).'], ...
    upper(func_name), var_name);
