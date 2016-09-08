function [new] = extractregion(old,polyname);
%   extractregion      - extracts a region in polygon from image using mouse or existing polygon
% [new] = extractregion(old,polyname);
%
% Function to extract a region from an image using the mouse 
% saves polygon mask  in 'tmp_lastregion'
% region can be inserted back with INSERTREGION
%
% use polyname if given as polygon
%

if nargin==0,help extractregion;return;end;

pixval off
fprintf(1,'Select a polygon using the mouse...\n');

   if nargin==2
      str=sprintf('load %s ', polyname);
      eval(str)
   else
      [BW,Xi,Yi] = roipoly;  %Xi and Yi are the polyn. vortices
   end
% convert binary image to indexed imag
   [BW,map]=gray2ind(BW);
% convert indexed image to intensity image
   BW = ind2gray(BW,map);
% create the mask
   BW(find(BW==0))=NaN;
% mask out the zero-values
   new=double(BW).*old;

   figure;
   imagesc(new) ; axis image
   save tmp_lastregion BW
