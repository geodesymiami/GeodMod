function [cmap,clim] = demcmap(varargin)
%DEMCMAP Colormaps appropriate to terrain elevation data
%
%  DEMCMAP(map) creates and assigns a colormap appropriate for elevation data.
%  The colormap has the number of land and sea colors in proportion to
%  the maximum elevations and depths in the matrix map. With no output
%  arguments, the colormap is applied to the current figure and the axes
%  CLim property is set so that the interface between the land and sea
%  is correct.
%
%  DEMCMAP(map,n) makes a colormap with a length of n.
%
%  DEMCMAP(map,n,rgbseamat,rgblandmat) allows the default colors for
%  sea and land to be replaced by colormap matrices. Default colors for either
%  land or sea can be retained by providing an empty matrix, [], in place
%  of a colormap matrix. The current figure colormap can be specified by
%  entering the string 'window' in place of either rgb matrix.  The input
%  colormap matrices may be of any size. The actual colors used in the
%  created colormap are interpolated.
%
%  DEMCMAP('mode',map,n,...) uses the 'mode' string to define the colormap
%  construction.  If 'mode' = 'size', n is the length of the colormap,
%  with a default value of n = 64.  If 'mode' = 'inc', n is the altitude
%  increment assigned to each color, with a default value of n = 100.  If omitted,
%  'mode' = 'size' is assumed.
%
%  [cmap,clim] = DEMCMAP(...) returns the colormap matrix and color axis vector,
%  but does not apply them to the current figure.

% Copyright 1996-2006 The MathWorks, Inc.
% $Revision: 1.13.4.4 $  $Date: 2006/12/10 20:04:59 $
% Written by:  A. Kim, W. Stumpf, E. Byrns

if nargin == 0
    errorOnMissingUI(mfilename)
end

if nargin == 1 && ~ischar(varargin{1})
    colorstr   = [];
    map        = varargin{1};
    sizearg    = [];
    rgbseamat  = [];
    rgblandmat = [];

elseif nargin == 2 && ~ischar(varargin{1})
    colorstr   = [];
    map        = varargin{1};
    sizearg    = real(varargin{2});
    rgbseamat  = [];
    rgblandmat = [];

elseif nargin == 4 && ~ischar(varargin{1})
    colorstr   = [];
    map        = varargin{1};
    sizearg    = real(varargin{2});
    rgbseamat  = varargin{3};
    rgblandmat = varargin{4};

elseif nargin == 2 && ischar(varargin{1})
    colorstr   = varargin{1};
    map        = varargin{2};
    sizearg    = [];
    rgbseamat  = [];
    rgblandmat = [];

elseif nargin == 3 && ischar(varargin{1})
    colorstr   = varargin{1};
    map        = varargin{2};
    sizearg    = real(varargin{3});
    rgbseamat  = [];
    rgblandmat = [];

elseif nargin == 5 && ischar(varargin{1})
    colorstr   = varargin{1};
    map        = varargin{2};
    sizearg    = real(varargin{3});
    rgbseamat  = varargin{4};
    rgblandmat = varargin{5};

else
    eid = sprintf('%s:%s:incorrectInputs', getcomp, mfilename);
    error(eid,'%s','Incorrect number of arguments');
end


%  Set the default color string
if isempty(colorstr)
    colorstr = 'size';
end

%  Input testing
if isempty(map)
    eid = sprintf('%s:%s:invalidMap', getcomp, mfilename);
    error(eid,'%s','Map input must not be empty.');
end

if ~isempty(sizearg)
    if length(sizearg) > 1
        eid = sprintf('%s:%s:', getcomp, mfilename);
        error(eid,'%s','N must be a scalar')
    elseif sizearg <= 0
        eid = sprintf('%s:%s:invalidN', getcomp, mfilename);
        error(eid,'%s','N must be positive')
    end
end

if ~isempty(rgbseamat) && ~strcmp(rgbseamat,'window')
    rgbseamat = ignoreComplex(rgbseamat, mfilename, 'RGBSEAMAT');
    if ndims(rgbseamat) > 2 || size(rgbseamat,2) ~= 3
        eid = sprintf('%s:%s:invalidColorMap', getcomp, mfilename);
        error(eid,'%s','Colormap inputs must be a 3 column matrix.')
    elseif any(rgbseamat(:) < 0) || any(rgbseamat(:) > 1)
        eid = sprintf('%s:%s:invalidColorMap', getcomp, mfilename);
        error(eid,'%s','Colormap inputs must be between 0 and 1')
    end
end

if ~isempty(rgblandmat) && ~strcmp(rgblandmat,'window')
    rgblandmat = ignoreComplex(rgblandmat, mfilename, 'RBGLANDMAT');
    if ndims(rgblandmat) > 2 || size(rgblandmat,2) ~= 3
        eid = sprintf('%s:%s:invalidColorMap', getcomp, mfilename);
        error(eid,'%s','Colormap inputs must be a 3 column matrix.')
    elseif any(rgblandmat(:) < 0) || any(rgblandmat(:) > 1)
        eid = sprintf('%s:%s:invalidColorMap', getcomp, mfilename);
        error(eid,'%s','Colormap inputs must be between 0 and 1')
    end
end

%  Compute the colormap
switch colorstr
    case 'size'
        [cmap0,clim0] = demcmap1(map,floor(sizearg),rgbseamat,rgblandmat);
    case 'inc'
        [cmap0,clim0] = demcmap2(map,sizearg,rgbseamat,rgblandmat);
    otherwise
        eid = sprintf('%s:%s:invalidParam', getcomp, mfilename);
        error(eid,'%s','Unrecognized color string')
end

%  Set the output arguments if necessary
if nargout==0;
    caxis(clim0);
    colormap(cmap0)
elseif nargout == 1;
    cmap = cmap0;
elseif nargout == 2;
    cmap = cmap0;
    clim = clim0;
end

%------------------------------------------------------------------
function [rgbmap,clim] = demcmap1(map,ncolors,rgbseamat,rgblandmat)

%DEMCMAP1  Generation of colormaps - fixed number of colors
%
%	Purpose
%		Creates colormaps for digital elevation maps (DEM's).  User can
%		provide number of desired colors and RGB color limit values for sea
%		and land, which can be any n-color (n-by-3) matrix.  A spectrum of
%		colors are built from the specified color values.  The default value
%		for ncolors is 64.  Default RGB color limits are provided, using
%		atlas-like colors.  User can specify the default 3-element sea or
%		3-element land color limits with the empty matrix [] in place of the
%		actual color values.
%
%	Synopsis
%
%		demcmap1(map)
%		demcmap1(map,ncolors)
%		demcmap1(map,ncolors,rgbseamat,rgblandmat)
%		[rgbmap,clim] = demcmap1(...)
%


if nargin == 1
    ncolors = [];
    rgbseamat = [];
    rgblandmat = [];
elseif nargin == 2
    rgbseamat = [];
    rgblandmat = [];
elseif nargin == 4
    if ~isempty(rgbseamat)
        hsvseamat = rgb2hsv(rgbseamat);
    end
    if ~isempty(rgblandmat)
        hsvlandmat = rgb2hsv(rgblandmat);
    end
else
    eid = sprintf('%s:%s:incorrectInputs', getcomp, mfilename);
    error(eid,'%s','Incorrect number of arguments')
end

%  Set defaults
if isempty(ncolors);   
    ncolors = 64;   
end

if isempty(rgbseamat)
    hsvseamat =	[2/3 1 0.2; 2/3 1 1; 0.5 1 1];				% dk b, b, c
elseif strcmp(rgbseamat,'window')
    hsvseamat = rgb2hsv(colormap);
end

if isempty(rgblandmat)
    hsvlandmat = [5/12 1 0.4; 0.25 0.2 1; 5/72 1 0.4];		% g-b, lt y-g, br
elseif strcmp(rgblandmat,'window')
    hsvlandmat = rgb2hsv(colormap);
end

%  Make sure at least two colors are selected
if ncolors <= 1;  
    ncolors = 2;  
end

minval = min(map(:));
maxval = max(map(:));
if minval == maxval;   
    maxval = minval+1;  
end

cmn = minval;
cmx = maxval;

% determine appropriate number of sea and land colors
if minval>=0
    nsea = 0;
    nland = ncolors;
elseif maxval<=0
    nland = 0;
    nsea = ncolors;
else
    % find optimal ratio of land to sea colors
    maxminratio = maxval/abs(minval);
    n1 = floor(ncolors/2);
    n2 = ceil(ncolors/2);
    if maxminratio>1
        sea = (1:n1)';
        land = (ncolors-1:-1:n2)';
    else
        land = (1:n1)';
        sea = (ncolors-1:-1:n2)';
    end
    ratio = land./sea;
    errors = abs(ratio - maxminratio) / maxminratio;
    indx = find(errors==min(min(errors)));
    nsea = sea(indx);
    nland = land(indx);

    % determine color limits
    seaint = abs(minval)/nsea;
    landint = maxval/nland;
    if seaint>=landint
        interval = seaint;
    else
        interval = landint;
    end
    cmn = -nsea*interval*(1 + 1e-9);		% zero values treated as land
    cmx = nland*interval;
end

% generate colormap for sea
nseamat = length(hsvseamat(:,1));
if nseamat==1 && nsea>1
    hsvseamat = [hsvseamat;hsvseamat];
    nseamat = 2;
end

if nsea==0
    hsvsea = [];
elseif nsea>0 && nsea<=nseamat
    temp = flipud(hsvseamat);
    hsvsea = flipud(temp(1:nsea,:));
else
    nseabands = nseamat - 1;
    nleft = nsea - nseamat;
    if nseabands==1
        nseavec = nleft;
    else
        for n=1:nseabands-1
            nseavec(n,1) = round(nleft/nseabands);
        end
        nseavec(nseabands,1) = nleft - sum(nseavec);
        nseavec = sort(nseavec);
    end
    hsvsea = hsvseamat(1,:);
    for n=1:nseabands
        deltahsv = (hsvseamat(n+1,:)-hsvseamat(n,:)) / (nseavec(n)+1);
        low = hsvseamat(n,:) + deltahsv;
        high = hsvseamat(n+1,:);
        clear hsvtemp
        for m=1:3
            if deltahsv(m)==0
                hsvtemp(:,m) = low(m)*ones(nseavec(n)+1,1);
            else
                hsvtemp(:,m) = (low(m):deltahsv(m):high(m))';
            end
        end
        hsvsea = [hsvsea; hsvtemp];
    end
end

% generate colormap for land
nlandmat = length(hsvlandmat(:,1));
if nlandmat==1 && nland>1
    hsvlandmat = [hsvlandmat;hsvlandmat];
    nlandmat = 2;
end
if nland==0
    hsvland = [];
elseif nland>0 && nland<=nlandmat
    hsvland = hsvlandmat(1:nland,:);
else
    nlandbands = nlandmat - 1;
    nleft = nland - nlandmat;
    if nlandbands==1
        nlandvec = nleft;
    else
        for n=1:nlandbands-1
            nlandvec(n,1) = round(nleft/nlandbands);
        end
        nlandvec(nlandbands,1) = nleft - sum(nlandvec);
        nlandvec = flipud(sort(nlandvec));
    end
    hsvland = hsvlandmat(1,:);
    for n=1:nlandbands
        deltahsv = (hsvlandmat(n+1,:)-hsvlandmat(n,:)) / (nlandvec(n)+1);
        low = hsvlandmat(n,:) + deltahsv;
        high = hsvlandmat(n+1,:);
        clear hsvtemp
        for m=1:3
            if deltahsv(m)==0
                hsvtemp(:,m) = low(m)*ones(nlandvec(n)+1,1);
            else
                hsvtemp(:,m) = (low(m):deltahsv(m):high(m))';
            end
        end
        hsvland = [hsvland; hsvtemp];
    end
end

rgbs = hsv2rgb([hsvsea; hsvland]);

if nargout==0;	
    caxis([cmn cmx]); 
    colormap(rgbs)
else
    rgbmap = rgbs;   
    clim = [cmn cmx];
end

%--------------------------------------------------------------------
function [rgbmap,clim] = demcmap2(map,surfdatainc,rgbseamat,rgblandmat)

%DEMCMAP2  Generation of colormaps - surface data increment
%
%	Purpose
%		Creates colormaps for digital elevation maps (DEM's).  User can
%		provide the surface data increment and RGB color limit values for sea
%		and land, which can be any n-color (n-by-3) matrix.  A spectrum of
%		colors are built from the specified color values.  The default value
%		for surfdatainc is 100.  Default RGB color limits are provided,  using
%		atlas-like colors.  User can specify the default 3-element sea or
%		3-element land color limits with the empty matrix [] in place of the
%		actual color values.
%
%	Synopsis
%
%		demcmap2(map)
%		demcmap2(map,surfdatainc)
%		demcmap2(map,surfdatainc,rgbseamat,rgblandmat)
%		[rgbmap,clim] = demcmap2(...)


if nargin == 1
    surfdatainc = [];     
    rgbseamat = [];    
    rgblandmat = [];
elseif nargin == 2
    rgbseamat = [];   
    rgblandmat = [];
elseif nargin == 4
    if ~isempty(rgbseamat)
        hsvseamat = rgb2hsv(rgbseamat);    
    end
    if ~isempty(rgblandmat)
        hsvlandmat = rgb2hsv(rgblandmat);  
    end
else
    eid = sprintf('%s:%s:incorrectInputs', getcomp, mfilename);
    error(eid,'%s','Incorrect number of arguments')
end

%  Set defaults
if isempty(surfdatainc)
    surfdatainc = 64;  
end
if isempty(rgbseamat)
    hsvseamat =	[2/3 1 0.2; 2/3 1 1; 0.5 1 1];				% dk b, b, c
elseif strcmp(rgbseamat,'window')
    hsvseamat = rgb2hsv(colormap);
end
if isempty(rgblandmat)
    hsvlandmat = [5/12 1 0.4; 0.25 0.2 1; 5/72 1 0.4];		% g-b, lt y-g, br
elseif strcmp(rgblandmat,'window')
    hsvlandmat = rgb2hsv(colormap);
end

minval = min(map(:));
maxval = max(map(:));
if minval == maxval;  
    maxval = minval+1; 
end

% determine appropriate number of sea and land colors
% determine color limits
if minval>=0
    nsea = 0;
    lowland = floor(minval/surfdatainc);
    highland = ceil(maxval/surfdatainc);
    nland = highland - lowland;

    cmn = lowland*surfdatainc;
    cmx = highland*surfdatainc;

elseif maxval<=0
    nland = 0;
    shallowsea = floor(abs(minval)/surfdatainc);
    deepsea = ceil(abs(minval)/surfdatainc);
    nsea = deepsea - shallowsea;

    cmn = -deepsea*surfdatainc;
    cmx = -shallowsea*surfdatainc;

else
    nsea = ceil(abs(minval)/surfdatainc);
    nland = ceil(maxval/surfdatainc);

    cmn = -nsea*surfdatainc*(1 + 1e-9);	% zero values treated as land;
    cmx = nland*surfdatainc;

end

% generate colormap for sea
nseamat = length(hsvseamat(:,1));
if nseamat==1
    hsvseamat = [hsvseamat;hsvseamat];
    nseamat = 2;
end
if nsea==0
    hsvsea = [];
elseif nsea>0 && nsea<=nseamat
    temp = flipud(hsvseamat);
    hsvsea = flipud(temp(1:nsea,:));
else
    nseabands = nseamat - 1;
    nleft = nsea - nseamat;
    if nseabands==1
        nseavec = nleft;
    else
        for n=1:nseabands-1
            nseavec(n,1) = round(nleft/nseabands);
        end
        nseavec(nseabands,1) = nleft - sum(nseavec);
        nseavec = sort(nseavec);
    end
    hsvsea = hsvseamat(1,:);
    for n=1:nseabands
        deltahsv = (hsvseamat(n+1,:)-hsvseamat(n,:)) / (nseavec(n)+1);
        low = hsvseamat(n,:) + deltahsv;
        high = hsvseamat(n+1,:);
        clear hsvtemp
        for m=1:3
            if deltahsv(m)==0
                hsvtemp(:,m) = low(m)*ones(nseavec(n)+1,1);
            else
                hsvtemp(:,m) = (low(m):deltahsv(m):high(m))';
            end
        end
        hsvsea = [hsvsea; hsvtemp];
    end
end

% generate colormap for land
nlandmat = length(hsvlandmat(:,1));
if nlandmat==1
    hsvlandmat = [hsvlandmat;hsvlandmat];
    nlandmat = 2;
end
if nland==0
    hsvland = [];
elseif nland>0 && nland<=nlandmat
    hsvland = hsvlandmat(1:nland,:);
else
    nlandbands = nlandmat - 1;
    nleft = nland - nlandmat;
    if nlandbands==1
        nlandvec = nleft;
    else
        for n=1:nlandbands-1
            nlandvec(n,1) = round(nleft/nlandbands);
        end
        nlandvec(nlandbands,1) = nleft - sum(nlandvec);
        nlandvec = flipud(sort(nlandvec));
    end
    hsvland = hsvlandmat(1,:);
    for n=1:nlandbands
        deltahsv = (hsvlandmat(n+1,:)-hsvlandmat(n,:)) / (nlandvec(n)+1);
        low = hsvlandmat(n,:) + deltahsv;
        high = hsvlandmat(n+1,:);
        clear hsvtemp
        for m=1:3
            if deltahsv(m)==0
                hsvtemp(:,m) = low(m)*ones(nlandvec(n)+1,1);
            else
                hsvtemp(:,m) = (low(m):deltahsv(m):high(m))';
            end
        end
        hsvland = [hsvland; hsvtemp];
    end
end

rgbs = hsv2rgb([hsvsea; hsvland]);

if nargout==0;	
    caxis([cmn cmx]); 
    colormap(rgbs)
else
    rgbmap = rgbs;    
    clim = [cmn cmx];
end
