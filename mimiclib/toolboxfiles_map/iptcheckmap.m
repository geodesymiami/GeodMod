function iptcheckmap(map, function_name, variable_name, argument_position)
%IPTCHECKMAP Check validity of colormap.
%   IPTCHECKMAP(MAP,FUNC_NAME,VAR_NAME,ARG_POS) checks to see if
%   MAP is a valid MATLAB colormap and issues a formatted error
%   message if it is invalid. 
%
%   FUNC_NAME is a string that specifies the name used in the formatted
%   error message to identify the function checking the colormap.
%
%   VAR_NAME is a string that specifies the name used in the formatted
%   error message to identify the argument being checked.
%
%   ARG_POS is a positive integer that indicates the position of
%   the argument being checked in the function argument list. 
%   IPTCHECKMAP converts this number to an ordinal number and includes
%   this information in the formatted error message.
%
%   Example
%   -------
%    
%       bad_map = ones(10);
%       iptcheckmap(bad_map,'func_name','var_name',2)
%
%   See also IPTCHECKHANDLE, IPTCHECKINPUT, IPTCHECKNARGIN, IPTCHECKSTRS,
%            IPTNUM2ORDINAL.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $  $Date: 2006/06/15 20:10:17 $

if ~isa(map,'double') || isempty(map) || ...
      (ndims(map) ~= 2) || (size(map,2) ~= 3) || ...
      issparse(map)
    msgId = sprintf('Images:%s:badMapMatrix',function_name);
    error(msgId,'Function %s expected its %s input argument, %s, to be a valid colormap.\nValid colormaps must be nonempty, double, 2-D matrices with 3 columns.', ...
          upper(function_name), iptnum2ordinal(argument_position), variable_name);
end

if any(map(:) < 0) || any(map(:) > 1)
    msgId = sprintf('Images:%s:badMapValues',function_name);
    error(msgId,'Function %s expected its %s input argument, %s, to be a valid colormap.\nValid colormaps cannot have values outside the range [0,1].', ...
          upper(function_name), iptnum2ordinal(argument_position), variable_name);
end
