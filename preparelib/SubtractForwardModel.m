% THIS SHOULD BE CHANGED TO BE WORKING AS A PROPER FUNCTION
  function [igram]=SubtractForwardModel(igram,opt)
%
% function of Data2Igram if opt.SubtractModel given in min file)
%
%                  - calculates LOS signal for a forward model using the satellite
%                    info from igram and one set of model parameters given in a file
%                    with same format as the usual geodmod output model.mat, 
%                    with the filename given in the *.min as 'opt.SubtractModel'. 
%                    The model parameters are assumed to be m/yr
%                  - Subtracts the forward model LOS signal from igram.data
%                    can be used e.g. to subtract interseismic signal from
%                    a variety of igrams (forward model signal gets adjusted 
%                    from m/yr to the actual timeperiod of the igram dates.
%       
%       INPUT:
%
%		igram(i)	  
%       opt        - providing name of 'opt.SubtractModel'
%
%
%       OUTPUT:
%
%       igram      - where igram.data is filtered, and igram.unfiltereddata is the
%                    original data (see end of script)
%     
%  Tini 08/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HAS TO BE ADJUSTED SO THAT IT WORKS FOR IGRAM(i).......... (see add_mask_to_igram.m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate radarlook angle from igram information:
for i=1:length(igram)

rearth       = 6.36e6; 
incanglevec  = asin((igram(i).sat_height+rearth)/rearth *sin(igram(i).los.look_ref4/180*pi)); 
radarlook    = [sin(incanglevec)*cos(igram(i).heading) -sin(incanglevec)*sin(igram(i).heading) cos(incanglevec)];

%--------------------------------------------------------------------------

% Calculate x_posting and y_posting for km grid spacing from igram 
% (as done for basemap in Dem2Basemap.m)

if strcmp(igram(i).x_unit,'degres') | strcmp(igram(i).x_unit,'degrees')
    x_last    = igram(i).x_first + igram(i).x_step*size(igram(i).data,2);
    y_last    = igram(i).y_first + igram(i).y_step*size(igram(i).data,1);

    x_center  = (igram(i).x_first+x_last)/2.;
    y_center  = (igram(i).y_first+y_last)/2.;

    pt1       = [y_center,igram(i).x_first] ;
    pt2       = [y_center,         x_last ] ;
    x_posting = distance(pt1,pt2,almanac('earth','WGS84')) / size(igram.data,2) ;

    pt1       = [igram(i).y_first,x_center] ;
    pt2       = [         y_last ,x_center] ;
    y_posting = distance(pt1,pt2,almanac('earth','WGS84'))/ size(igram.data,1);
else
    error('user error: x_unit %s not recognized -- exiting', igram.x_unit) ;
end

% Make 2 x n coordinate array for use in ForwardModel 

  x_length = size(igram(i).data,2)*x_posting;
  y_length = size(igram(i).data,1)*y_posting;
  x        = linspace(0,x_length,size(igram(i).data,2));
  y        = linspace(0,y_length,size(igram(i).data,1));
  [xx,yy]  = meshgrid(x,y);
  coord    = [xx(:) yy(:)]';
  
%--------------------------------------------------------------------------

% Get Model parameters for ForwardModel (format as stored in model.mat)
% (Change to read in directly from template file??)

S =load(opt.SubtractModelFile);
modelopt = S.modelopt;

%--------------------------------------------------------------------------

% Run ForwardModel to get interseismic displacement field in LOS direction
[m,n] = size(coord);
radarlook = repmat(radarlook(:),1,n);
[iu_radarlook,iu,iu_sum] = ForwardModel(modelopt.par.xy,coord,radarlook,modelopt);

%--------------------------------------------------------------------------

% Reshape vector into igram array format
% I hope the additions to igram changing the length of this structure 
% dont mess up anything else in the code....

[m,n]         = size(igram(i).data);
interseismic  = reshape(iu_radarlook,[m,n]);

%--------------------------------------------------------------------------

% Change the rate in the interseismic model from m/yr to time-period of
% igram:

interseismic  = interseismic * ((igram(i).t2-igram(i).t1)/365); % double

%--------------------------------------------------------------------------

% Subtract the model from the igram

igram(i).data            = igram(i).data - interseismic;

% save(out_name,'igram'); % not necessary since done in Data2Igram.m
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

end
disp(sprintf('SubtractForwardModel is working jippie'))
return % to Data2Igram



















