function [rgbim,clim] = PlotRGB_Demjpg_N(colim,dem,shadefact,colm,ampsc,shadsat,colsat,colorbaropt,CLim);
%    PlotRGB_Demjpg_N    - Creates a color image with shaded topo background.
%
% usage:  [rgbim,clim] = PlotRGB_Demjpg_N(colim,dem,shadefact,colm,ampsc,shadsat,colsat,file);
%
% Input:
% colim          - (n,m) data matrix, give the color in output image
% dem            - (n,m,3) shading (see preparedem.pl)
% shadefact  - (1,1) Scale the ration data/shade (default = 1)
% colm           - (63,3) colormap matrix, default is 'jet'
% ampsc          - (1,1) to scale the brightness (default=1), try 1.5!
% shadsat        - (2,1) vector with shading saturation values,
%                    default: max,min in shading (e.g. try [-50 50])
% colsat         - (2,1) vector with color saturation values,
%                    default: max,min of input color matrix
% colorbaropt    - option structure for colorbar
%                  colorbaropt.ColorbarParam [xcenter ycenter width height]
%                  colorbaropt.ColorbarTitle 
%                  colorbaropt.ColorbarLabels [# of labels]
%
% Output:
% rgbim          - (n,m,3) rgb image
% clim           - (2,1) vector with minimum and maximum colorvalues
%
% Sigurjon Jonsson, 7 July 2004
% Modify Noel Gourmelen, September 2005
% Falk Amelung Oct 31, allos
%
% Set default values, depending on number of input arguments
%if nargin<5;  sun_direction=[20 45];  end
if nargin<3;  shadefact=1;  end
if nargin<4;  colm=jet;  end
if nargin<5;  ampsc=1;   end
if nargin<6;  shadsat=1; end
if nargin<7;  colsat=1;  end
if nargin<8; colorbaropt=[];  end
%colsat=[20;-20];
%------------------------------------------------------------
% Check dem:

if isfloat(dem)==0 ; dem=double(dem) ;  end
if max(dem(:))>1 ;   dem=dem/max(dem(:)) ;  end
    
    
% find number of lines and column in data-matrix
[ll,cc] = size(colim);

%----------------------------------------------------------------------
% Prepare color image data

% Saturate color values, or find max,min of colorimage
if length(colsat)>1
  minc = colsat(1);
  maxc = colsat(2);
  colim(colim<minc)=minc;
  colim(colim>maxc)=maxc;
else
  minc = min(colim(:));
  maxc = max(colim(:));
end

% Report min,max elevation values
clim = [minc maxc]';


% sbaker: need this to work for contour plots
%  clim=CLim;
%  minc=CLim(1);
%  maxc=CLim(2);





% adjust dem values, values between 0 and 63, as
% we'll keep top-bin empty for coloring the NaNs gray
d = colim-minc;
%d = ( d/max(d(:)) )*63;
d  = ( d/(maxc-minc) )*63;
d = round(d);            % has to be integers

% Here we modify the colortable, as the NaNs are colored
% by the top value, put gray for that
kjet = [colm(1:size(colm,1)-1,:); 0.9 0.9 0.9];

% Transform DEM image to rgb image using 'colm' colormap
rgbim = ind2rgb(d,kjet);

% Scale with amplitudes of RGB image using the scaled shading values
%whos dem
%shadefact=0.4
%TODO: Here we should work with uint8. For that we should have the shaded relief as
%      uint8 carried through the code. add_shade2Igram can be easily modified
%      but the code here needs to be modified as well so that the multiplication works

demfact=dem(:,:,1).^shadefact ;
for k=1:3
  out  = rgbim(:,:,k);   % re-use 'out'
  out  = out.*(demfact);
  %rgbim(:,:,k) = out.*dem(:,:,k);    % FA Jan 2008: re-instated old version (next line). Don't know why this was put in.
  rgbim(:,:,k) = out;
end


nan_pt = find(isnan(colim(:)));

% replacement RGB value at NaN point to Reflectance value
   R_sc_at_nan_pt = (nan_pt);
   [nrows ncols] = size(colim);

   clear colim R_sc_at_nan_pt amp

   rgbim(nan_pt)                 = dem(nan_pt)             ; % R value
   rgbim(nan_pt+(nrows*ncols)*1) = dem(nan_pt+(nrows*ncols)*1) ; % G value
   rgbim(nan_pt+(nrows*ncols)*2) = dem(nan_pt+(nrows*ncols)*2) ; % B value

   rgbim = uint8(ceil(rgbim*255));  %FA 10/2008

% convert to unsigned integer to save memory
% TODO: Conversions could be done at an earlier stage to further save memory
%   rgbim = uint8(ceil(rgbim*255));
   
%imagesc(bounds(1):bounds(2):bounds(3),bounds(4):bounds(5):bounds(6),rgbim);hold on;
%axis xy





