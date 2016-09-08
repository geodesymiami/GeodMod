function [igram_name,timeseries_name,motion_name,qt_name,dataset_name,suffix]=GenerateNames(dir_in,dir_out,motion2qtopt)
%GenerateQtorGridName   - generates name for quadtree or grid file
%
%usage: [igram_name,motion_name,qt_name,dataset_name,suffix]=GenerateNames(dir_in,dir_out,motion2qtopt)
%
%Input:   motion_name     name to which suffic will be appended  (e.g. motion_RsatA3 --> motion_RsatA3_qt)
%         opt           option structure containing opt.suffix
%
%Output:  outname       new name       (e.g. motion_RsatA3_qt)
%         suffix        generated according to the following rule:
%                       if opt not given      'Grid'
%                       if opt.suffix given   'opt.suffix'
%                       else use name according to opt.Method
%         
%  Falk Amelung, September 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 [junk,junk,satbeam]=extract_ProjectName(dir_in);

 igram_name      = [ dir_out filesep 'igram_'      junk satbeam '.mat'] ;
 timeseries_name = [ dir_out filesep 'timeseries_' junk satbeam '.mat'] ;
 motion_name     = [ dir_out filesep 'motion_'     junk satbeam '.mat'] ;

  if ~exist('motion2qtopt','var')
     suffix='Grid';
  else
      if isfield(motion2qtopt,'suffix') 
         suffix=motion2qtopt.suffix;
      else  
         if isfield(motion2qtopt,'Method')
            if     strmatch('Quadtree',motion2qtopt.Method)        suffix='Qt'   ;      
            elseif strmatch('Grid',motion2qtopt.Method)            suffix='Grid' ;     
            elseif strmatch('ProjectProfile',motion2qtopt.Method)  suffix='Profile' ;   
            end                          
         else
            suffix='';
         end
      end
   end
  
  [pathstr, name, ext] = fileparts(motion_name); 
  qt_name=[pathstr filesep  name '_' suffix ext];
  dataset_name=[pathstr filesep 'dataset_' suffix ext];


