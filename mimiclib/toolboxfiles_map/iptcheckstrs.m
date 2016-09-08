function out = iptcheckstrs(in, valid_strings, function_name, ...
                         variable_name, argument_position)
%IPTCHECKSTRS Check validity of text string.
%   OUT = IPTCHECKSTRS(IN,VALID_STRINGS,FUNC_NAME,VAR_NAME,ARG_POS) checks 
%   the validity of the text string IN. If the text string matches one of
%   the text strings in the cell array VALID_STRINGS, IPTCHECKSTRS returns
%   the valid text string in OUT. If the text string does not match, 
%   IPTCHECKSTRS issues a formatted error message.
%
%   IPTCHECKSTRS looks for a case-insensitive nonambiguous match between
%   IN and the strings in VALID_STRINGS.
%
%   VALID_STRINGS is a cell array containing text strings.
%
%   FUNC_NAME is a string that specifies the name used in the formatted
%   error message to identify the function checking text strings.
%
%   VAR_NAME is a string that specifies the name used in the formatted
%   error message to identify the argument being checked.
%
%   ARG_POS is a positive integer that indicates the position of
%   the argument being checked in the function argument list. 
%   IPTCHECKSTRS converts this number to an ordinal number and includes
%   this information in the formatted error message.
%
%   Example
%   -------
%       % To trigger this error message, define a cell array of some text
%       % strings and pass in another string that isn't in the cell array. 
%       iptcheckstrs('option3',{'option1','option2'},'func_name','var_name',2)
%
%   See also IPTCHECKHANDLE, IPTCHECKINPUT, IPTCHECKMAP, IPTCHECKNARGIN
%            IPTNUM2ORDINAL.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2006/10/04 22:39:00 $

% Except for IN, input arguments are not checked for validity.

% If you do not wish to specify argument position in resulting error
% messages, provide [] for ARG_POS.

iptcheckinput(in, {'char'}, {'row'}, function_name, variable_name, ...
              argument_position);

matches = strncmpi(in,valid_strings,numel(in));
num_matches = sum(matches);

if num_matches == 1
  out = valid_strings{matches};

else
  out = substringMatch(valid_strings(matches));
    
  if isempty(out)
    % Convert valid_strings to a single string containing a space-separated list
    % of valid strings.
    list = '';
    for k = 1:length(valid_strings)
      list = [list ', ' valid_strings{k}];
    end
    list(1:2) = [];
    
    if ~isempty(argument_position)
      msg1 = sprintf('Function %s expected its %s input argument, %s,', ...
                     upper(function_name), iptnum2ordinal(argument_position), ...
                     variable_name);
    else
      msg1 = sprintf('Function %s expected its input argument %s', ...
                     upper(function_name), variable_name);
    end
    
    msg2 = 'to match one of these strings:';
    
    if num_matches == 0
      msg3 = sprintf('The input, ''%s'', did not match any of ', in);
      msg4 = 'the valid strings.';
      id = sprintf('Images:%s:unrecognizedStringChoice', function_name);
      
    else
      msg3 = sprintf('The input, ''%s'', matched more than one ', in);
      msg4 = 'valid string.';
      id = sprintf('Images:%s:ambiguousStringChoice', function_name);
    end
    
    error(id,'%s\n%s\n\n  %s\n\n%s%s', msg1, msg2, list, msg3,msg4);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = substringMatch(strings)
%   STR = substringMatch(STRINGS) looks at STRINGS (a cell array of
%   strings) to see whether the shortest string is a proper substring of
%   all the other strings.  If it is, then substringMatch returns the
%   shortest string; otherwise, it returns the empty string.

if isempty(strings)
  str = '';
else
  len = cellfun('prodofsize',strings);
  [tmp,sortIdx] = sort(len);
  strings = strings(sortIdx);
  
  start = regexpi(strings(2:end), ['^' strings{1}]);
  if isempty(start) || (iscell(start) && any(cellfun('isempty',start)))
    str = '';
  else
    str = strings{1};
  end
end
