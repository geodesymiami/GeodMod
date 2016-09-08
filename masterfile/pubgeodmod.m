function []=pubgeodmod(fname,outdir)
%pubgeodmod    - generates an html file based on a *min file 
%
%  usage:  pubgeodmod file.min [dir] 
%
%    e.g.: geodmod Lengai_1fault.min  
%          pubgeodmod Lengai_1fault.min '/RAID1/amelung/models/html'
%
%          Options given 
%
%          'dir_in'            Input directory (contains e.g. RsatA3)       [default 'off']
%          'plotdataopt'       plot options                                 [default 'off' ]
%                        
%  Saves data as igram structure.
%  If corfile is given it is used 
%  
%  V1.0  Falk Amelung, October 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [pathstr,name,ext] = fileparts(fname) ;
 fname_for_publish = fullfile(pathstr,[name '.m']) ;
 fname_html        = fullfile(pathstr,[name '.html']) ;

 if exist('outdir','var') 
        opts.outputDir = outdir; 
        full_fname_html = fullfile(opts.outputDir,fname_html) ;
    else 
        opts.outputDir = name; 
        full_fname_html = fullfile(pwd,opts.outputDir,fname_html) ;
 end

 disp('To view html:')
 disp(['firefox ' full_fname_html]); 

 if ~exist(opts.outputDir,'dir')  mkdir(opts.outputDir) ; end
 %logmessage(sprintf('output directory: %s ',opts.outputDir));

 fid = fopen(fname_for_publish,'w');
     fprintf(fid,'type %s \n',fname) ;
     fprintf(fid,'geodmod %s \n',fname) ;
 fclose(fid) ;
 
 disp('publishing to html....')
 disp('when finished, view with:')
 disp(['firefox ' full_fname_html]); 
 publish(fname_for_publish,opts);

 delete(fname_for_publish)
