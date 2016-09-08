function n = ecc2n(parm)
%ECC2N  n-value of ellipse with given eccentricity
%
%  n = ECC2N(mat) computes the parameter n of an ellipse (or
%  ellipsoid of revolution) given the eccentricity.  n is defined
%  as (a-b)/(a+b), where a = semimajor axis, b = semiminor axis.
%  If the input is a column vector, then each column element is assumed
%  to be an eccentricity.  If the input has two columns, then the
%  second column is assumed to be the eccentricity.  This allows
%  geoid vectors from ALMANAC to be used as inputs.  If the
%  input is a n x m matrix, where m ~= 2, then each element is assumed
%  to be an eccentricity and the corresponding n is calculated.
%
%  See also N2ECC, ECC2FLAT, MAJAXIS, MINAXIS.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.10.4.4 $  $Date: 2007/11/09 20:23:30 $
% Written by:  E. Byrns, E. Brown

%  Dimension tests
if min(size(parm)) == 1 && ndims(parm) <= 2
	col = min(size(parm,2), 2);   % Select first or second column
	eccent = parm(:,col);         % First col if single vector input
	                              % Second col if two column inputs (eg. geoid vecs)
else
    eccent = parm;                %  General matrix input
end

%  Ensure real inputs
eccent = ignoreComplex(eccent, mfilename, 'eccentricity');

%  Compute n
n = (1 - sqrt(1 - eccent.^2)) ./ (1 + sqrt(1 - eccent.^2));
