 function [stdMap] = makeSTD(data,win)

%
% Compute the std of a input matrice over the area defined by the window
% sixe "win"
%
% 	[stdMap] = makeSTD(data,win)
% 
% data     :      input structure. (LoadData format) 
%
% win is pixel radius of circle window
%
%
% N. Gourmelen Jan. 2010
%

%%%%%%%%%%%%%%%%%%%%
if ~isstruct(data) data.data = data ;  end

se1   = strel('disk',win) ;
nhood = getnhood(se1)       ;

stdMap = stdfilt(data.data,nhood) ;


