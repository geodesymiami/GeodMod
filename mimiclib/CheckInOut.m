function [out_flag]=CheckInOut(in_list,out_list);
%CheckInOut      -  checks for existence of input and output file 
%
%usage: [out_flag]=CheckInOut(in_file,out_file);
%
%Input:   in_file, outfile    
%
%Output:  out_flag  (true if file exists, else false)
%
%  skips checking if string is empty (e.g. in_file='')
%
%  Falk Amelung, September 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[st,i]=dbstack ;
%logmessage(sprintf('',2));

%
% First check out_list
%
if ~iscell(out_list)  out_list=cellstr(out_list);  end
   out_flag=false;
   for i=1:length(out_list) 
       out_file=out_list{i};
       if ~isempty(out_file)
          if  exist(out_file,'file') ||  exist([out_file '.mat'],'file') 
              out_flag=true;
              logmessage(sprintf('out_file found: %s --- skipping %s',out_file,st(2).name),2);
              return
          else
              out_flag=false;
              logmessage(sprintf('will generate: %s',out_file),2);
          end
       end
   end

%
% Now check infiles 
%
if ~iscell(in_list)  in_list=cellstr(in_list);  end
   for i=1:length(in_list) 
       in_file=in_list{i};
       if ~isempty(in_file)
          if  exist(in_file,'file') 
              logmessage(sprintf('in_file  found: %s',in_file),2);
          else
              logmessage(sprintf('exiting --- file %s does not exist',in_file));
              error(sprintf('user error %s->%s: exiting --- file %s does not exist',st(2).name,mfilename,in_file));
          end
      end
end

end
