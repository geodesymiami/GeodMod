function   [OutStruct]=add_struct(struct1,struct2)

%
%   Function to add one structure (struct2) to another one (struct1)
%
%       [OutStruct]=add_struct(struct1,struct2)
%
% N Gourmelen - August 2007
%

if nargin < 2
    exit ('Error, not enough arguments');
end

fields=fieldnames(struct2);

for ni=1:length(fields)
    
    struct1.(genvarname(list2str(fields(ni))))=struct2.(genvarname(list2str(fields(ni))));
    
end

OutStruct=struct1;
