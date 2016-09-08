function []=logplot(plot_command,out_name,varargin);
%logplot                - generates plot and saves it to a file 
%
%usage: []=logplot('plot_igram(igram)');
%
%  function plots the plotstring and writes it to a file
%  if outname is given
%  Future versions will write to plot file depending on the level. A threshold level 
%  will be given as a global variable in the same way as for logmessage
%
%  if plot_commmand is empty the current figure will be saved. 
%  this is for compatibility with PlotModel and can be removed once PlotModel is obsolete
%
%  Falk Amelung, September 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
   error ('logplot needs at least 3 arguments --- exiting')
end

global plot_visible plot_save 

  if ~isempty(plot_command)                    % do not open a figure and plot if plot_command is an empty string

     figure;
     if ~plot_visible
      set(gcf,'visible','off');    %Not yet testes
     end

     if strcmp(plot_command,'patch')        % ToDo: call patch from a higher-level program such as PlotPatch that generates the figure handle, axis
  %      figure
  %      if ~plot_visible
  %         set(gcf,'visible','off');    %Not yet testes
  %      end
        feval(plot_command,varargin{:}); 
        axis ij ; axis image;               % June 2007 change. Was axis xy.
        
%      elseif strcmp(plot_command,'PlotData')
%         %feval({'[rgbim,kml_string] =' plot_command},varargin{:});
%         %feval(plot_command,varargin{:});
%         [rgbim]=feval(plot_command,varargin{:});
%         
%         
%         % save kmz file if rgbim and kml_string is returned
%         if exist(kml_string,'var')
% 			imwrite(rgbim,'image_overlay.jpg','JPEG')
% 			kml_file = 'doc.kml';
% 			image_file = 'image_overlay.jpg';
% 			filename = strcat(out_name,'.kmz');
% 		% Write the contents of kml_string to kml_file
% 			fid = fopen( kml_file, 'wt'); 
% 			fprintf(fid, kml_string );
% 			fclose(fid);
% 		% Create the KMZ File 
% 			zip(filename,{kml_file,image_file});  % zip the kml file and image jpeg
% 			zipped_file = strcat(filename,'.zip'); % the zip command appends .zip to the file name which we don't want
% 			movefile(zipped_file,filename); % move the zipped file to the name ending with .kmz
% 			delete(kml_file,image_file); % delete the kml and image file to save space since they are in the .kmz file now
%         end
     
     else
 
        feval(plot_command,varargin{:}); 

     end

  end

%
% return if out_name is empty
%
  if isempty(out_name) return ; end

%
% save plotfile
%
    if strcmp(plot_save,'pdf')
      [pathstr, name, ext] = fileparts(out_name);
      plotdir              = fullfile (pathstr,'plot');  
      nout_name            = fullfile (plotdir,[name '.' plot_save]);
      if ~exist(plotdir,'dir')   mkdir(plotdir);  end
      
      logmessage (sprintf('generating plot: %s',nout_name))   
      
      tmp=dbstack;
      surface_flag=sum(strmatch('PlotSurface3D',{tmp.name}));
      
      if  surface_flag                               % there is an issue with transparency rendering. It does not work out of terminal emulation mode (matlabnn)
           %print('-dpdf','-painter',nout_name)
      else  
           print('-dpdf',nout_name);
      end
  end
