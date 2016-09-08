function iptchecknargin(low, high, numInputs, function_name)
%IPTCHECKNARGIN Check number of input arguments.
%   IPTCHECKNARGIN(LOW,HIGH,NUM_INPUTS,FUNC_NAME) checks whether
%   the number of input arguments NUM_INPUTS is in the range specified
%   by LOW and HIGH. If NUM_INPUTS is not in this range, IPTCHECKNARGIN
%   issues a formatted error message.
%
%   LOW must be a scalar nonnegative integer.
%
%   HIGH must be a scalar nonnegative integer or Inf.
%
%   FUNC_NAME is a string that specifies the name used in the formatted
%   error message to identify the function checking its input
%   arguments.
%
%   Example
%   -------
%
%   Create a function and use IPTCHECKNARGIN to check that the 
%   number of arguments passed to the function is within the 
%   expected range.
%
%       function test_function(varargin)
%       iptchecknargin(1,3,nargin,mfilename);
%
%   Trigger the error message by executing the function at 
%   the MATLAB command line, specifying more the expected 
%   number of arguments.
%   
%       test_function(1,2,3,4)
%  
%   See also IPTCHECKHANDLE, IPTCHECKINPUT, IPTCHECKMAP, IPTCHECKSTRS,
%            IPTNUM2ORDINAL.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2004/11/19 17:36:50 $

% Input arguments are not checked for validity.

if numInputs < low
  msgId = sprintf('Images:%s:tooFewInputs', function_name);
  if low == 1
    msg1 = sprintf('Function %s expected at least 1 input argument', ...
                   upper(function_name));
  else
    msg1 = sprintf('Function %s expected at least %d input arguments', ...
                   upper(function_name), low);
  end
  
  if numInputs == 1
    msg2 = 'but was called instead with 1 input argument.';
  else
    msg2 = sprintf('but was called instead with %d input arguments.', ...
                   numInputs);
  end
  
  error(msgId, '%s\n%s', msg1, msg2);
  
elseif numInputs > high
  msgId = sprintf('Images:%s:tooManyInputs', function_name);

  if high == 1
    msg1 = sprintf('Function %s expected at most 1 input argument', ...
                   upper(function_name));
  else
    msg1 = sprintf('Function %s expected at most %d input arguments', ...
                   upper(function_name), high);
  end
  
  if numInputs == 1
    msg2 = 'but was called instead with 1 input argument.';
  else
    msg2 = sprintf('but was called instead with %d input arguments.', ...
                   numInputs);
  end
  
  error(msgId, '%s\n%s', msg1, msg2);
end

  
    