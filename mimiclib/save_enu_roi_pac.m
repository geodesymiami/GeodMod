function []=save_enu_roi_pac(filename,data)
% save_enu_roi_pac    - saves a *enu file (EastNorthUp) with a rsc file
%
% function []=save_enu_roi_pac(filename,data)
%
% writes data in typical structure format (3 datasets) into enu file 
% similar as hgt format but 3 interleaved datasets
% NaNs are replaced by 0s.
%
% TODO: Could be simply modified to write a *hgt file if the
% if length(data)=2, and also to take the format ('float32', 'int16') as
% input option with 'float32' as the default.
%
if nargin==0,help save_enu_roi_pac;return;end;

% filenames
enufname=[filename '.enu'];
rscfname=[filename '.enu' '.rsc'];

%
% save data 
%
logmessage(sprintf('saving... %s', enufname))
fid1=fopen(enufname,'wb');
   data(1).data(isnan(data(1).data))=0;
   data(2).data(isnan(data(2).data))=0;
   data(3).data(isnan(data(3).data))=0;

   Ttmp=[data(1).data';data(2).data';data(3).data'];
   fwrite(fid1,Ttmp,'float32');
fclose(fid1);

%
% save rscfile 
%
str=char(     ['WIDTH                 ' num2str(size(data(1).data,2))]) ;    
str=char(str, ['FILE_LENGTH           ' num2str(size(data(1).data,1))]) ;    
str=char(str, ['X_FIRST               ' num2str(data(1).x_first)]) ;    
str=char(str, ['Y_FIRST               ' num2str(data(1).y_first)]) ;    
str=char(str, ['X_STEP                ' num2str(data(1).x_step)]) ;    
str=char(str, ['Y_STEP                ' num2str(data(1).y_step)]) ;    
str=char(str, ['X_UNIT                ' num2str(data(1).x_unit)]) ;    
str=char(str, ['UNIT                  ' num2str(data(1).Unit)]) ;    

logmessage(sprintf('saving... %s', rscfname))
fid1=fopen(rscfname,'w');
for i=1:size(str,1)              
    fprintf(fid1,'%s\n',str(i,:));
end
fclose(fid1);
