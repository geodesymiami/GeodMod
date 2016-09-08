function  [cam_position]=viewdirections3d_generate(modelopt) 
%
%   viewdirections3d_generate  - generates a cell string with view directions for PlotSurface3D
%

cam_position      = {[sind(20+180)*100 cosd(20+180)*100 10]}; %default
   
if modelopt.N_disloc || modelopt.N_multidisloc
   strike         =   modelopt.par.xy(5);
   cam_position   = {[sind(strike+180)*100 cosd(strike+180)*100 -modelopt.par.xy(3)]   [sind(strike+90)*100 cosd(strike+90)*100 -modelopt.par.xy(3)] cam_position{:}} ;
end
