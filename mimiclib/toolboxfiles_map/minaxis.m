function semiminor = minaxis(in1,in2)
%MINAXIS  Semiminor axis of ellipse with given semimajor axis and
%         eccentricity
%
%  b = MINAXIS(semimajor,e) computes the semiminor axis of an ellipse
%  (or ellipsoid of revolution) given the semimajor axis and eccentricity.
%  The input data can be scalar or matrices of equal dimensions.
%
%  b = MINAXIS(vec) assumes a 2 element vector (vec) is supplied,
%  where vec = [semimajor, e].
%
%  See also AXES2ECC, FLAT2ECC, MAJAXIS, N2ECC.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.4 $  $Date: 2007/11/09 20:24:46 $
% Written by:  E. Byrns, E. Brown

error(nargchk(1, 2, nargin, 'struct'))

if nargin == 1
    if ~isequal(sort(size(in1)),[1 2])
        error(['map:' mfilename ':mapError'], ...
            'Input must be a 2 element vector')
    else
        in1 = ignoreComplex(in1, mfilename, 'vec');
        semimajor = in1(1);
        eccent    = in1(2);
    end
elseif nargin == 2
    if ~isequal(size(in1),size(in2))
        error(['map:' mfilename ':mapError'], ...
            'Inconsistent input dimensions')
    else
        semimajor = ignoreComplex(in1, mfilename, 'semimajor');
        eccent    = ignoreComplex(in2, mfilename, 'eccentricity');
    end
end

%  Compute the semiminor axis
semiminor = semimajor .* sqrt(1 - eccent.^2);
