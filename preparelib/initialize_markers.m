function [marker1D,marker2D]  = initialize_markers(marker1D,marker2D,basemap)
%initialize_multidislocopt   - converts multidisloc option string into structure
 
  if iscell(marker1D)
    for i=1:length(marker1D)
       if ~isfield(marker1D{i},'xy') && isfield(marker1D{i},'lola')   
           marker1D{i}.xy        =lola2xy(str2num(marker1D{i}.lola{1})',basemap,1);  
           marker1D{i}.symb      = marker1D{i}.lola{2};
           marker1D{i}.txt       = marker1D{i}.lola{3};
           marker1D{i}.lola        = str2num(marker1D{i}.lola{1})';
       end
    end
  else
    marker1D=false;
  end

  
    if iscell(marker2D)
    for i=1:length(marker2D)
       if ~isfield(marker2D{i},'xy') && isfield(marker2D{i},'lola')   
           marker2D{i}.xy        =lola2xy(str2num(marker2D{i}.lola{1})',basemap,1);  
           marker2D{i}.symb      = marker2D{i}.lola{2};
           marker2D{i}.txt       = marker2D{i}.lola{3};
           marker2D{i}.lola      = str2num(marker2D{i}.lola{1})';
       end
    end
else
    marker2D=false;
end
 

