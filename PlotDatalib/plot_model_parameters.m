function [s] = plot_model_parameters(modelopt,x_unit,varargin);
%   plot_mode_parameters    - plots dislocations and magma sources
% usage:  [] = plot_model_parameters(model,x_unit,PlotType);
%
% FA, May 2007, 
% FA, Oct 2008, adopted for 2D,3D and 1D (CrossSection)
% TODO:  there is lots of unnecessary legacy code. Streamlining necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     if length(varargin)==1
         PlotType = varargin{1};
     else
         Profile  = varargin{1};
         PlotType = varargin{2};
     end
     
     if strcmp(PlotType,'1D') && length(varargin)==1 PlotType='2D'; end     % needed because too much plotdata code is run in PlotData in any case
     
     if strcmp(PlotType,'3D') load slipcol; end    % Load colormap

     N_disloc         = modelopt.N_disloc;
     N_mogi           = modelopt.N_mogi;
     N_penny          = modelopt.N_penny;
     N_mctigue        = modelopt.N_mctigue;
     N_pCDM           = modelopt.N_pCDM;
     N_yang           = modelopt.N_yang;
     N_multidisloc    = modelopt.N_multidisloc;
     N_visco1d        = modelopt.N_visco1d;
     N_lockedandcreep = modelopt.N_lockedandcreep;

     if strcmp(x_unit,'degres') | strcmp(x_unit,'degrees')
        par = modelopt.par.lola;
        if strcmp(PlotType,'3D') errordlg('3D plotting currently only for x_unit=km, not degrees'); error(' user error -- exiting'); end
     else
        par = modelopt.par.xy;
     end

     if N_multidisloc  
         N_disloc = modelopt.multidislocopt.N_disloc;
         par      = multidislocpar2dislocpar(par,modelopt.multidislocopt,x_unit);
     end

     
     if N_disloc>=1
          switch PlotType
          case {'2D' 'Map'}
               if N_disloc >= 1   plot(par(6),par(7),'k*');  dislocplot(par(1:10),x_unit); dislocplot_surfaceintersection(par(1:10),x_unit,PlotType); par(1:10)=[]; end
               if N_disloc >= 2   plot(par(6),par(7),'k*');  dislocplot(par(1:10),x_unit); dislocplot_surfaceintersection(par(1:10),x_unit,PlotType); par(1:10)=[]; end  % FA 3/2010: need a control flag for plotting of
               if N_disloc >= 3   plot(par(6),par(7),'k*');  dislocplot(par(1:10),x_unit); dislocplot_surfaceintersection(par(1:10),x_unit,PlotType); par(1:10)=[]; end  %  surface projection (depends on HorzOff of multidislocopt)
               if N_disloc >= 4   plot(par(6),par(7),'k*');  dislocplot(par(1:10),x_unit); dislocplot_surfaceintersection(par(1:10),x_unit,PlotType); par(1:10)=[]; end
          case {'3D'}
               pm              = reshape(par(1:N_disloc*10)',10,N_disloc)';
               i               = find(pm(:,10)~=0);           % opening dislocations
               ind             = zeros(size(pm));
               ind(i,:)        = 1;                           % all indices of opening dislcoations=1
               i_op            = find(ind==1);
               i_mag           = find(ind~=1);
               [s_op]          = plot_disloc3d(par(i_op),'op');
               [s_mag]         = plot_disloc3d(par(i_mag),'mag');
               s               = [0 max([max(s_op) max(s_mag)] )];   % used for colorbar only
               par(1:N_disloc*10) = [];
          case {'1D' 'CrossSection'}
                
                dislocplot_surfaceintersection(par(1:10),x_unit,Profile,PlotType)
                
                for i=1:N_disloc                 
                  [fx,fy,fz]          = flakes(par(1:7));
                  disloc_profile_dist = Profile.unitvec*[fx'-Profile.xy(1,1); fy'-Profile.xy(1,2)]; 
                  plot(disloc_profile_dist,-fz,'-r','LineWidth',5);
                  par(1:10)           =[];
               end
          end
     end
     
     if N_mogi>=1
          switch PlotType
          case {'2D' 'Map'}
               if N_mogi       >= 1  plot(par(1),par(2),'k+');  par(1:4) =[]; end
               if N_mogi       >= 2  plot(par(1),par(2),'k+');  par(1:4) =[]; end
               if N_mogi       >= 3  plot(par(1),par(2),'k+');  par(1:4) =[]; end
          case {'3D'}
               par_mogi        = par(1:N_mogi*4);
               [s]             = plot_mogi3d(par_mogi);
               s               = [0 abs(s)];
               par(1:N_mogi*4) = [] ;
          case {'1D' 'CrossSection'}  
               logmessage('--> Need to write CrossSection plot code for this source')
          end
     end
     
     if N_penny>=1 
          switch PlotType
          case {'2D' 'Map'}
               if N_penny      >= 1  plot(par(1),par(2),'k+');  par(1:5) =[]; end
               if N_penny      >= 2  plot(par(1),par(2),'k+');  par(1:5) =[]; end
          case {'3D'}                 
                par_penny = par(1:N_penny*5);
                [s]       = plot_penny3d(par_penny);
                s         = [0 abs(s)];
                par(1:N_penny*5)  = [];
          case {'1D' 'CrossSection'}
                logmessage('--> Need to write CrossSection plot code for this source')
          end
     end

     if N_mctigue>=1 
          switch PlotType
          case {'2D' 'Map'}
               if N_mctigue    >= 1  plot(par(1),par(2),'k+');  par(1:5) =[]; end
               if N_mctigue    >= 2  plot(par(1),par(2),'k+');  par(1:5) =[]; end  
          case {'3D'}
               logmessage('--> Need to write 3D plot program for this source')
               s=0.1;
          case {'1D' 'CrossSection'}
               logmessage('--> Need to write CrossSection plot code for this source')
          end
     end

     % copied from Mogi
     if N_pCDM>=1
          switch PlotType
          case {'2D' 'Map'}
               if N_mogi       >= 1  plot(par(1),par(2),'k+');  par(1:4) =[]; end
               if N_mogi       >= 2  plot(par(1),par(2),'k+');  par(1:4) =[]; end
               if N_mogi       >= 3  plot(par(1),par(2),'k+');  par(1:4) =[]; end
          case {'3D'}
               par_mogi        = par(1:N_mogi*4);
               [s]             = plot_mogi3d(par_mogi);
               s               = [0 abs(s)];
               par(1:N_mogi*4) = [] ;
          case {'1D' 'CrossSection'}  
               logmessage('--> Need to write CrossSection plot code for this source')
          end
     end
     
     if N_yang>=1
          switch PlotType
          case {'2D' 'Map'}
               if N_yang       >= 1  plot(par(1),par(2),'k+');  par(1:5) =[]; end
               if N_yang       >= 2  plot(par(1),par(2),'k+');  par(1:5) =[]; end  
          case {'3D'}
               logmessage('--> Need to write 3D plot program for this source')
               s=0.1;
          case {'1D' 'CrossSection'}
               logmessage('--> Need to write CrossSection plot code for this source')
          end
     end       
         
             
     
   

    
