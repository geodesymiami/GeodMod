function []=logmessage(string,level);
%logmessage      -  displays message and writes to log file  
%
%usage: []=logmessage('Hello friend');
%
%  This function is meant to record messages to
%  various  log files depending on the message level.
%  The writing of log files will be switched 'on'/'off'
%  by a global variable. A time stamp could also be controlled
%  by a global variable.
%
%  In the current version the string is only written to the screen.
%
%  Falk Amelung, September 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[st,i]=dbstack ;
if length(st) >= 3 
   str=sprintf('%-13s>%-13s: %s', st(3).name,st(2).name,string);
elseif length(st)==2
   str=sprintf('%-26s: %s',st(2).name,string);
else
   str=sprintf('%-26s:', mfilename,string);
end
disp(str)
