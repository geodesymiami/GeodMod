%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script to find all dependencies of FUNCTIONSTR
% copies all entires matching matchdir1,2,3 
% into the destination directory allfunctions
% FIRST RUN SETBASEDIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
functionstr={'masterfile/MakeDataModel' }

matchdir1  = '/home/amelung/matlab'    ;
matchdir2  = '/RAID6/insar_lab/matlab' ;
matchdir3  = basedir ;
%matchdir3  = '/qqq' ;      comment out if basedir entries not desired

destdir  = [basedir '/allfunctions/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:length(functionstr)
   flist=depfun(functionstr{j});            %TODO: There must be a simple way to exclude all built-in functions
   flist{1}=[];
   ind=zeros(length(flist),1);

   % keep functions containing certain strings (e.g. 'insar_lab','amelung')
      for i=1:length(flist)
         if isempty(strfind(flist{i},matchdir1)) & isempty(strfind(flist{i},matchdir2)) & isempty(strfind(flist{i},matchdir3))   ind(i) = 1; end
      end
      flist(find(ind==1))=[];

   % remove functions that match certain string, e.g. basedir, 'insar_lab'
      %ind=zeros(length(flist),1);
      %for i=1:length(flist)
      %   if ~isempty(strfind(char(flist(i)),basedir)) ind(i) = 1; end
      %end
      %flist(find(ind==1))=[];

   % generate destination list by replacing matlabhomedir by destinationdir in the file path
   for i=1:length(flist)
        if ~isempty(strfind(flist{i},matchdir1))  tmpchar=strrep(flist{i},matchdir1,destdir);  end
        if ~isempty(strfind(flist{i},matchdir2))  tmpchar=strrep(flist{i},matchdir2,destdir);  end
        if ~isempty(strfind(flist{i},matchdir3))  tmpchar=strrep(flist{i},matchdir3,destdir);  end
        [pathstr,name,ext,versn] = fileparts(tmpchar) ;
        dlist{i}=pathstr;
   end
   disp(flist)

   % copy files into destdir
   for i=1:length(flist)
       mkdir(dlist{i})            ;
       copyfile(flist{i},dlist{i});
   end
end

