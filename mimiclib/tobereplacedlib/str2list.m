function [strlist]=str2list(str)
% str2list     -  converts  string into a cellstr
%
% usage:   strlist=str2list(str)
%
% strlist contains elements of str separated by blanks 
%
% FA May 2002
%
%  see also:   list2str
%
str=strjust(str,'left'); str=deblank(str);   %removes heading and trailing blanks
str(findstr(str,'  '))=[];
list=[1 findstr(str,' ')+1];

for i=1:length(list)
    if i<length(list)
       tmp=str(list(i):list(i+1)-2);
    else
       tmp=str(list(i):end);
    end
    %strlist{i}=tmp;          % changed jun20 fa
    strlist{i}=deblank(tmp);
end




