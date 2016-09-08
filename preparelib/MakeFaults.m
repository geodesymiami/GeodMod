%Script to convert GMT fault files into *mat file as for ML2002-2005.min
% First replace all section delimiters by
% NaN  NaN
% then do the following:

% THIS CODE SHOULD GO INTO PREPARE SO THAT WE CAN READ *txt FILES. THE SEPARATOR COULD BE GIVEN AS AN OPTION

fid = fopen('/RAID6/insarlab/geophys_data/Haiti/Lines/caribbean_faults.lola.txt');
      C     = textscan(fid,'%f %f');
      long  = [C{1}'];
      lat   = [C{2}'];
      Lllh  = [long' lat'];
save('/RAID6/insarlab/geophys_data/Haiti/Lines/caribbean_faults.mat');
