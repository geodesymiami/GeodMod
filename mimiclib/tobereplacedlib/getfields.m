function [outfield] = getfields(structure,field,outtype);

%
% [outfield] = getfields(structure,field);
%
%   structure : input structure
%
%   field     : field for which you want the values
%
%   outtype   : 'cell', 'string' (default)
%
%
% N. Gourmelen, November 2007
%

for ni=1:length(structure)  
    tmp=structure(ni).(genvarname(field)); 
    if (nargin==3 && strcmp(outtype,'cell'))
        tmp={tmp};
    end
    outfield(ni,:)=tmp;
end

    