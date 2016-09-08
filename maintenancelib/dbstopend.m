function []=dbstopend(fname);
%dbstopend      -  invokes debugger at last line of file  
%
%usage: []=dbstopend('PreparePlot');
%
%  Falk Amelung, June 2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fullname=which(fname);
[q,str]=unix(['wc -l ' fullname]) ;
cstr=str2list(str);
str=['dbstop in ' fname ' at ' cstr{1}];
disp(str)
eval(str)
