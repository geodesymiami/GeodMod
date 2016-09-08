function [out_flag]=CheckToolboxes();
%CheckToolboxes      -  checks for existence of matlab toolboxes required by geodmod
%
%usage: [out_flag]=CheckToolboxes();
%
%Output:  out_flag  (true if all toolboxes  exist,  false otherwise)
%
%  Falk Amelung, October 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list = {'image_toolbox' 'map_toolbox' 'optimization_toolbox' 'satistics_toolbox'};

out_flag= true;

for i=1:length(list);
    if ~license('test',list{1})
       txt_str=['geodmod needs: ' list{i} '. You may see errors ....'] ;
       errordlg(txt_str); error([txt_str 'user error -- exiting']); 
       out_flag=false;
       return
   end
end 
