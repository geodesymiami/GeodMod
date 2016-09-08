function [str]=list2str(cellstr)
% list2str     -  converts  cellstr into str (elements separated by blanks)
%
% usage:   str=list2str(cellstr)
%
% str contains elements of cellstr separated by blanks 
%
% FA May 2002
%
%  see also:   str2list
%
tmp=char(cellstr);
str=[];
for i=1:size(tmp,1)
  tmpstr=sprintf('%s?',cellstr{i}) ; str=strcat(str,tmpstr);
end

str=strrep(str,'?',' ');
