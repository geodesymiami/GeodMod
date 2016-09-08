function [outdata] = smooth(indata,indistance,smoothing_distance)

%
% Smooth data using a window averaging set by the keyword smoothing distance
%
% indata             : Data to smooth
% indistance         : Location of data
% smoothing_distance : window length (same unit as indistance)
%
% Output:
%            outdata: Smoothed data, same size as indata
%
% N. Gourmelen, Feb 2008
%

ld         = length(indata)       ;
outdata    = zeros(ld,1)          ;   
filterhalf = smoothing_distance/2 ; 

for ni=1:ld,
    position         = indistance(ni)                                    ;
    neighborid       = find(abs(position-indistance) <= filterhalf)      ;
    scale            = filterhalf - abs(position-indistance(neighborid))     ;
    outdata(ni)       = sum(indata(neighborid).*scale)./sum(scale)    ;
end
