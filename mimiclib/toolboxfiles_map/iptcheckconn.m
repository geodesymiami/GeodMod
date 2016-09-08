function iptcheckconn(conn,function_name,variable_name,arg_position)
%IPTCHECKCONN Check validity of connectivity argument.
%   IPTCHECKCONN(CONN,FUNC_NAME,VAR_NAME,ARG_POS) checks if CONN
%   is a valid connectivity argument. If it is invalid, the function
%   issues a formatted error message.
%
%   A connectivity argument can be one of the following scalar 
%   values, 1, 4, 6, 8, 18, or 26. A connectivity argument can also
%   be a 3-by-3-by- ... -by-3 array of 0s and 1s. The central element
%   of a connectivity array must be nonzero and the array must be
%   symmetric about its center. 
%
%   FUNC_NAME is a string that specifies the name used in the formatted
%   error message to identify the function checking the connectivity
%   argument.
%
%   VAR_NAME is a string that specifies the name used in the formatted
%   error message to identify the argument being checked.
%
%   ARG_POS is a positive integer that indicates the position of
%   the argument being checked in the function argument list. 
%   IPTCHECKCONN converts this number to an ordinal number and includes
%   this information in the formatted error message.
%
%   Class Support
%   -------------
%   CONN must be of class double or logical and must be real and nonsparse.
%
%   Example
%   -------   
%       % Create a 4-by-4 array and pass it as connectivity argument.
%       iptcheckconn(eye(4), 'func_name','var_name',2)  
%
%   See also IPTNUM2ORDINAL.

%    Copyright 1993-2004 The MathWorks, Inc.
%    $Revision: 1.1.8.1 $  $Date: 2004/08/10 01:50:45 $

iptcheckinput(conn,{'double' 'logical'},{'real' 'nonsparse'},...
              function_name,variable_name,arg_position);

if all(size(conn) == 1)
    if (conn ~= 1) && (conn ~= 4) && (conn ~= 8) && (conn ~= 6) && ...
                (conn ~= 18) && (conn ~= 26)

        msgId = sprintf('Images:%s:badScalarConn', function_name);
        msg1 = first_line(variable_name, function_name, arg_position);
        msg2 = 'A scalar connectivity specifier must be 1, 4, 6, 8, 18, or 26.';
        error(msgId,'%s\n%s',msg1,msg2);
    end
else
    if any(size(conn) ~= 3)
        msgId = sprintf('Images:%s:badConnSize', function_name);
        msg1 = first_line(variable_name, function_name, arg_position);
        msg2 = 'A nonscalar connectivity specifier must be 3-by-3-by- ... -by-3.';
        error(msgId,'%s\n%s',msg1,msg2);
    end
    
    if any((conn(:) ~= 1) & (conn(:) ~= 0))
        msgId = sprintf('Images:%s:badConnValue', function_name);
        msg1 = first_line(variable_name, function_name, arg_position);
        msg2 = 'A nonscalar connectivity specifier must contain only 0s and 1s.';
        error(msgId,'%s\n%s',msg1,msg2);
    end
    
    if conn((end+1)/2) == 0
        msgId = sprintf('Images:%s:badConnCenter', function_name);
        msg1 = first_line(variable_name, function_name, arg_position);
        msg2 = 'The central element of a connectivity specifier must be nonzero.';
        error(msgId,'%s\n%s',msg1,msg2);
    end
    
    if ~isequal(conn(1:end), conn(end:-1:1))
        msgId = sprintf('Images:%s:nonsymmetricConn', function_name);
        msg1 = first_line(variable_name, function_name, arg_position);
        msg2 = 'A connectivity specifier must be symmetric about its center.';
        error(msgId,'%s\n%s',msg1,msg2);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = first_line(variable_name, function_name, arg_position)

str = sprintf('Function %s expected its %s input argument, %s,\nto be a valid connectivity specifier.', ...
              upper(function_name), iptnum2ordinal(arg_position), variable_name);
