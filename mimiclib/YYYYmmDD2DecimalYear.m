function [DecYear]=YYYYmmDD2DecimalYear(datestring);
% YYYYmmDD2DecimalYear  - Converts a datestring in roi_pac format to decimal year (e.g. 
%
%  Falk Amelung, December 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    for i=1:size(datestring,1)
        year       = str2num(datestring(i,1:4))    ;
        month      = str2num(datestring(i,5:6))    ;
        day        = str2num(datestring(i,7:8))    ;
        DecYear(i) = year+(datenum(year,month,day) - datenum(year,1,1) ) / 365.25;
    end

    DecYear = DecYear';
