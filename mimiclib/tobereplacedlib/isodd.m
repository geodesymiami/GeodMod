function res = isodd(in)
% ISODD  True for odd numbers.
%   ISODD(IN) returns 1 if all elements of IN are odd,
%   and 0 otherwise.
%   See also ISEVEN, ISINT, REM, MOD.
%

% $Revision: 1.1.1.1 $  $Date: 2008-08-12 16:22:13 $
% Bert Kampes, 31-Mar-2000

res=0;
if(rem(in,2))res=1;end;

%%% EOF
