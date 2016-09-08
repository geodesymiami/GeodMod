function [list]=bracketedstring2cell(instring)
%bracketedstring2cell  - splits a brackedted string separated by blanks into a cell  
%
%usage: [list]=bracketed_dates_2_two_dates(instring);
%
%Input:   instring      string of the form '[ 030101 080115 ]'
%
%Output:  date1,date2   strings of the form '030101' '080115'
%
%  Falk Amelung, January 2009
%
%  TODO: It actually should return dates in the form [2003.12 2008.18] as this is the form used for time_span
%        and then we have to make SelectPairs working with this.
%  TODO: Even better would be a function that reads any input related to a date and converts in a format
%        such as 2003.12 used throughout geodomod.  Even better, we should use throughout geodmod the same matlab date
%        format.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  tmp           = strtrim(instring(strfind(instring,'[')+1:strfind(instring,']')-1));
  list          = str2list(tmp);
