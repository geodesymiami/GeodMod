function flat = ecc2flat(parm)
%ECC2FLAT  Flattening of ellipse with given eccentricity
%
%  f = ECC2FLAT(mat) computes the flattening of an ellipse (or
%  ellipsoid of revolution) given the eccentricity.  If the
%  input is a column vector, then each input is assumed to be an
%  eccentricity.  If the input has two columns, then the second
%  column is assumed to be the eccentricity.  This allows geoid
%  vectors from ALMANAC to be used as inputs.  If the input
%  is an n x m matrix, where m,n > 1, then each element is assumed
%  to be an eccentricity.
%
%  See also FLAT2ECC, ECC2N, MAJAXIS, MINAXIS.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.4 $  $Date: 2007/11/09 20:23:29 $
% Written by:  E. Byrns, E. Brown

%  Dimension tests
if min(size(parm)) == 1 && ndims(parm) <= 2
	col = min(size(parm,2), 2);   % Select first or second column
	eccent = parm(:,col);         % First col if single vector input
	                              % Second col if two column inputs (eg. geoid vecs)
else
    eccent = parm;        %  General matrix input
end

%  Ensure real inputs
eccent = ignoreComplex(eccent, mfilename, 'eccentricity');

%  Compute the flattening
flat = 1 - sqrt(1 - eccent.^2);
