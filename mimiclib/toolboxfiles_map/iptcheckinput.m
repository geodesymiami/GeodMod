%IPTCHECKINPUT Check validity of array.
%   IPTCHECKINPUT(A,CLASSES,ATTRIBUTES,FUNC_NAME,VAR_NAME, ARG_POS) checks
%   the validity of the array A and issues a formatted error message if it
%   is invalid. A can be an array of any class.
%
%   CLASSES is a cell array of strings containing the set of classes that A
%   is expected to belong to.  For example, if you specify CLASSES as
%   {'logical','cell'}, A is required to be either a logical array or a cell
%   array. The string 'numeric' is interpreted as an abbreviation for the
%   classes uint8, uint16, uint32, int8, int16, int32, single, double.
%
%   ATTRIBUTES is a cell array of strings specifying the set of attributes
%   that A must satisfy.  For example, if you specify ATTRIBUTES as {'real',
%   'nonempty','finite'}, A must be real and nonempty, and it must contain
%   only finite values.  The supported list of attributes includes:
%   
%       2d               nonempty            nonzero       row 
%       column           nonnan              odd           scalar
%       even             nonnegative         positive      twod
%       finite           nonsparse           real          vector       
%       integer        
%
%   FUNC_NAME is a string that specifies the name used in the formatted
%   error message to identify the function checking the input.
%
%   VAR_NAME is a string that specifies the name used in the formatted
%   error message to identify the argument being checked.
%
%   ARG_POS is a positive integer that indicates the position of
%   the argument being checked in the function argument list. 
%   IPTCHECKINPUT converts this number to an ordinal number and includes
%   this information in the formatted error message.
%
%   Example
%   -------
%       % To trigger this error message, create a three dimensional array
%       % and then check for the attribute '2d'.
%       A = [ 1 2 3; 4 5 6 ];
%       B = [ 7 8 9; 10 11 12];
%       C = cat(3,A,B);
%       iptcheckinput(C,{'numeric'},{'2d'},'my_func','my_var',2)
%
%   See also IPTCHECKHANDLE, IPTCHECKMAP, IPTCHECKNARGIN, IPTCHECKSTRS, 
%            IPTNUM2ORDINAL.

%   Copyright 1993-2007 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $  $Date: 2007/05/10 13:47:12 $

% If you do not wish to specify argument position in resulting error
% messages, provide [] for ARG_POS.