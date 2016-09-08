function outstats = regionprops(varargin)
%REGIONPROPS Measure properties of image regions (blob analysis).
%   STATS = REGIONPROPS(L,PROPERTIES) measures a set of properties for each
%   labeled region in the label matrix L. Positive integer elements of L
%   correspond to different regions. For example, the set of elements of L
%   equal to 1 corresponds to region 1; the set of elements of L equal to 2
%   corresponds to region 2; and so on. STATS is a structure array of
%   length max(L(:)). The fields of the structure array denote different
%   properties for each region, as specified by PROPERTIES. L can be
%   multidimensional.
%
%   STATS = REGIONPROPS(L,I,PROPERTIES) measures a set of properties for
%   each labeled region in the 2-D or N-D grayscale image I. L is a label
%   matrix that identifies the regions in I and must have the same
%   size as I. 
%
%   PROPERTIES can be a comma-separated list of strings, a cell array
%   containing strings, the string 'all', or the string 'basic'. The set of
%   valid measurement strings includes:
%
%   Shape measurements
%
%     'Area'              'EulerNumber'       'Orientation'               
%     'BoundingBox'       'Extent'            'Perimeter'          
%     'Centroid'          'Extrema'           'PixelIdxList' 
%     'ConvexArea'        'FilledArea'        'PixelList'
%     'ConvexHull'        'FilledImage'       'Solidity' 
%     'ConvexImage'       'Image'             'SubarrayIdx'            
%     'Eccentricity'      'MajorAxisLength' 
%     'EquivDiameter'     'MinorAxisLength'                   
%     
%   Pixel value measurements (requires grayscale image as an input)
%
%     'MaxIntensity'
%     'MeanIntensity'
%     'MinIntensity'
%     'PixelValues'
%     'WeightedCentroid'
%
%   Property strings are case insensitive and can be abbreviated.
%
%   If PROPERTIES is the string 'all', REGIONPROPS returns all of the
%   Shape measurements. If called with a grayscale image, REGIONPROPS also
%   returns Pixel value measurements. If PROPERTIES is not specified or if
%   it is the string 'basic', these measurements are computed: 'Area',
%   'Centroid', and 'BoundingBox'.
%
%   Perimeter should be used on a label matrix with contiguous regions, such
%   as L = bwlabel(BW). Otherwise, 'Perimeter' gives unexpected results on
%   discontiguous regions.
%
%   Note - REGIONPROPS and binary images
%   ------------------------------------
%   REGIONPROPS does not accept a binary image as its first input.  There
%   are two common ways to convert a binary image to a label matrix:
%
%       1.  L = bwlabel(BW);
%
%       2.  L = double(BW);
%
%   Suppose that BW were a logical matrix containing these values:
%
%       1 1 0 0 0 0
%       1 1 0 0 0 0
%       0 0 0 0 0 0
%       0 0 0 0 1 1
%       0 0 0 0 1 1
%
%   The first method of forming a label matrix, L = bwlabel(BW), results
%   in a label matrix containing two contiguous regions labeled by the
%   integer values 1 and 2.  The second method of forming a label matrix,
%   L = double(BW), results in a label matrix containing one
%   discontiguous region labeled by the integer value 1.  Since each
%   result is legitimately desirable in certain situations, REGIONPROPS
%   does not accept binary images and convert them using either method.
%   You should convert a binary image to a label matrix using one of
%   these methods (or another method if appropriate) before calling
%   REGIONPROPS.
%
%   Example
%   -------
%   Label the connected pixel components in the text.png image, compute
%   their centroids, and superimpose the centroid locations on the
%   image.
%
%       bw = imread('text.png');
%       L = bwlabel(bw);
%       s  = regionprops(L, 'centroid');
%       centroids = cat(1, s.Centroid);
%       imshow(bw)
%       hold on
%       plot(centroids(:,1), centroids(:,2), 'b*')
%       hold off
%
%   Class Support
%   -------------
%   The input label matrix L can have any numeric class and any dimension. L
%   must be real, nonsparse, and contain integers.
%
%   See also BWLABEL, BWLABELN, ISMEMBER, WATERSHED.

%   Copyright 1993-2007 The MathWorks, Inc.
%   $Revision: 1.7.4.13.4.1 $  $Date: 2007/12/28 21:38:57 $

[L,I,requestedStats,officialStats] = ParseInputs(varargin{:});

if ndims(L) > 2
    % Remove stats that aren't supported for N-D input and issue
    % warning messages as appropriate.
    requestedStats = PreprocessRequestedStats(requestedStats);
end

if isempty(requestedStats)
    eid = sprintf('Images:%s:noPropertiesWereSelected',mfilename);
    error(eid, 'No properties were selected.');
end

if (isempty(L))
    numObjs = 0;
else
    numObjs = round(double(max(L(:))));
end

% Initialize the stats structure array.
tempStats = {'PerimeterCornerPixelList'};
allStats = [officialStats; tempStats];
numStats = length(allStats);
empties = cell(numStats, numObjs);
stats = cell2struct(empties, allStats, 1);
% Initialize the statsAlreadyComputed structure array. Need to avoid
% multiple calculatations of the same property for performance reasons.
zz = cell(numStats, 1);
for k = 1:numStats
    zz{k} = 0;
end
statsAlreadyComputed = cell2struct(zz, allStats, 1);

% Calculate PixelIdxList
[stats, statsAlreadyComputed] = ...
    ComputePixelIdxList(L, stats, statsAlreadyComputed, numObjs);

% Compute statistics.
numRequestedStats = length(requestedStats);
for k = 1 : numRequestedStats
    switch requestedStats{k}

        case 'Area'
            [stats, statsAlreadyComputed] = ...
                ComputeArea(stats, statsAlreadyComputed);

        case 'FilledImage'
            [stats, statsAlreadyComputed] = ...
                ComputeFilledImage(L,stats,statsAlreadyComputed);

        case 'FilledArea'
            [stats, statsAlreadyComputed] = ...
                ComputeFilledArea(L,stats,statsAlreadyComputed);

        case 'ConvexArea'
            [stats, statsAlreadyComputed] = ...
                ComputeConvexArea(L, stats,statsAlreadyComputed);

        case 'Centroid'
            [stats, statsAlreadyComputed] = ...
                ComputeCentroid(L, stats, statsAlreadyComputed);

        case 'EulerNumber'
            [stats, statsAlreadyComputed] = ...
                ComputeEulerNumber(L,stats,statsAlreadyComputed);

        case 'EquivDiameter'
            [stats, statsAlreadyComputed] = ...
                ComputeEquivDiameter(L, stats, statsAlreadyComputed);

        case 'Extrema'
            [stats, statsAlreadyComputed] = ...
                ComputeExtrema(L,stats,statsAlreadyComputed);

        case 'BoundingBox'
            [stats, statsAlreadyComputed] = ...
                ComputeBoundingBox(L,stats,statsAlreadyComputed);

        case 'SubarrayIdx'
            [stats, statsAlreadyComputed] = ...
                ComputeSubarrayIdx(L,stats,statsAlreadyComputed);

        case {'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Eccentricity'}
            [stats, statsAlreadyComputed] = ...
                ComputeEllipseParams(L,stats,statsAlreadyComputed);

        case 'Solidity'
            [stats, statsAlreadyComputed] = ...
                ComputeSolidity(L,stats,statsAlreadyComputed);

        case 'Extent'
            [stats, statsAlreadyComputed] = ...
                ComputeExtent(L,stats,statsAlreadyComputed);

        case 'ConvexImage'
            [stats, statsAlreadyComputed] = ...
                ComputeConvexImage(L,stats,statsAlreadyComputed);

        case 'ConvexHull'
            [stats, statsAlreadyComputed] = ...
                ComputeConvexHull(L,stats,statsAlreadyComputed);

        case 'Image'
            [stats, statsAlreadyComputed] = ...
                ComputeImage(L,stats,statsAlreadyComputed);

        case 'PixelList'
            [stats, statsAlreadyComputed] = ...
                ComputePixelList(L,stats,statsAlreadyComputed);

        case 'Perimeter'
            [stats, statsAlreadyComputed] = ...
                ComputePerimeter(L,stats,statsAlreadyComputed);
            
        case 'PixelValues'
            [stats, statsAlreadyComputed] = ...
                ComputePixelValues(I,stats,statsAlreadyComputed);

        case 'WeightedCentroid'
            [stats, statsAlreadyComputed] = ...
                ComputeWeightedCentroid(L,I,stats,statsAlreadyComputed);

        case 'MeanIntensity'
            [stats, statsAlreadyComputed] = ...
                ComputeMeanIntensity(I,stats,statsAlreadyComputed);
            
        case 'MinIntensity'
            [stats, statsAlreadyComputed] = ...
                ComputeMinIntensity(I,stats,statsAlreadyComputed);
            
        case 'MaxIntensity'
            [stats, statsAlreadyComputed] = ...
                ComputeMaxIntensity(I,stats,statsAlreadyComputed);
    end
end

% Initialize the output stats structure array.
empties = cell(numRequestedStats, numObjs);
outstats = cell2struct(empties, requestedStats, 1); %#ok<NASGU>

fnames = fieldnames(stats);
deleteStats = fnames(~ismember(fnames,requestedStats));
outstats = rmfield(stats,deleteStats);

%%%
%%% ComputePixelIdxList
%%%
function [stats, statsAlreadyComputed] = ...
    ComputePixelIdxList(L, stats,statsAlreadyComputed,numobj)
%   A P-by-1 matrix, where P is the number of pixels belonging to
%   the region.  Each element contains the linear index of the
%   corresponding pixel.

statsAlreadyComputed.PixelIdxList = 1;

if ~isempty(L) && numobj ~= 0
    idxList = regionpropsmex(L, numobj);
    [stats.PixelIdxList] = deal(idxList{:});
end

%%%
%%% ComputeArea
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeArea(stats, statsAlreadyComputed)
%   The area is defined to be the number of pixels belonging to
%   the region.

if ~statsAlreadyComputed.Area
    statsAlreadyComputed.Area = 1;

    for k = 1:length(stats)
        stats(k).Area = size(stats(k).PixelIdxList, 1);
    end
end

%%%
%%% ComputeEquivDiameter
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeEquivDiameter(L, stats, statsAlreadyComputed)
%   Computes the diameter of the circle that has the same area as
%   the region.
%   Ref: Russ, The Image Processing Handbook, 2nd ed, 1994, page
%   511.

if ~statsAlreadyComputed.EquivDiameter
    statsAlreadyComputed.EquivDiameter = 1;

    if ndims(L) > 2
        NoNDSupport('EquivDiameter');
        return
    end

    [stats, statsAlreadyComputed] = ...
        ComputeArea(stats,statsAlreadyComputed);

    factor = 2/sqrt(pi);
    for k = 1:length(stats)
        stats(k).EquivDiameter = factor * sqrt(stats(k).Area);
    end
end

%%%
%%% ComputeFilledImage
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeFilledImage(L,stats,statsAlreadyComputed)
%   Uses imfill to fill holes in the region.

if ~statsAlreadyComputed.FilledImage
    statsAlreadyComputed.FilledImage = 1;

    [stats, statsAlreadyComputed] = ...
        ComputeImage(L,stats,statsAlreadyComputed);

    conn = conndef(ndims(L),'minimal');

    for k = 1:length(stats)
        stats(k).FilledImage = imfill(stats(k).Image,conn,'holes');
    end
end

%%%
%%% ComputeConvexArea
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeConvexArea(L,stats,statsAlreadyComputed)
%   Computes the number of "on" pixels in ConvexImage.

if ~statsAlreadyComputed.ConvexArea
    statsAlreadyComputed.ConvexArea = 1;

    if ndims(L) > 2
        NoNDSupport('ConvexArea');
        return
    end

    [stats, statsAlreadyComputed] = ...
        ComputeConvexImage(L,stats,statsAlreadyComputed);

    for k = 1:length(stats)
        stats(k).ConvexArea = sum(stats(k).ConvexImage(:));
    end
end

%%%
%%% ComputeFilledArea
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeFilledArea(L,stats,statsAlreadyComputed)
%   Computes the number of "on" pixels in FilledImage.

if ~statsAlreadyComputed.FilledArea
    statsAlreadyComputed.FilledArea = 1;

    [stats, statsAlreadyComputed] = ...
        ComputeFilledImage(L,stats,statsAlreadyComputed);

    for k = 1:length(stats)
        stats(k).FilledArea = sum(stats(k).FilledImage(:));
    end
end

%%%
%%% ComputeConvexImage
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeConvexImage(L,stats,statsAlreadyComputed)
%   Uses ROIPOLY to fill in the convex hull.

if ~statsAlreadyComputed.ConvexImage
    statsAlreadyComputed.ConvexImage = 1;

    if ndims(L) > 2
        NoNDSupport('ConvexImage');
        return
    end

    [stats, statsAlreadyComputed] = ...
        ComputeConvexHull(L,stats,statsAlreadyComputed);
    [stats, statsAlreadyComputed] = ...
        ComputeBoundingBox(L,stats,statsAlreadyComputed);

    for k = 1:length(stats)
        M = stats(k).BoundingBox(4);
        N = stats(k).BoundingBox(3);
        hull = stats(k).ConvexHull;
        if (isempty(hull))
            stats(k).ConvexImage = false(M,N);
        else
            firstRow = stats(k).BoundingBox(2) + 0.5;
            firstCol = stats(k).BoundingBox(1) + 0.5;
            r = hull(:,2) - firstRow + 1;
            c = hull(:,1) - firstCol + 1;
            stats(k).ConvexImage = roipoly(M, N, c, r);
        end
    end
end

%%%
%%% ComputeCentroid
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeCentroid(L, stats, statsAlreadyComputed)
%   [mean(r) mean(c)]

if ~statsAlreadyComputed.Centroid
    statsAlreadyComputed.Centroid = 1;

    [stats, statsAlreadyComputed] = ...
        ComputePixelList(L,stats,statsAlreadyComputed);

    for k = 1:length(stats)
        stats(k).Centroid = mean(stats(k).PixelList,1);
    end

end

%%%
%%% ComputeEulerNumber
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeEulerNumber(L,stats,statsAlreadyComputed)
%   Calls BWEULER on 'Image' using 8-connectivity

if ~statsAlreadyComputed.EulerNumber
    statsAlreadyComputed.EulerNumber = 1;

    if ndims(L) > 2
        NoNDSupport('EulerNumber');
        return
    end

    [stats, statsAlreadyComputed] = ...
        ComputeImage(L,stats,statsAlreadyComputed);

    for k = 1:length(stats)
        stats(k).EulerNumber = bweuler(stats(k).Image,8);
    end
end

%%%
%%% ComputeExtrema
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeExtrema(L, stats, statsAlreadyComputed)
%   A 8-by-2 array; each row contains the x and y spatial
%   coordinates for these extrema:  leftmost-top, rightmost-top,
%   topmost-right, bottommost-right, rightmost-bottom, leftmost-bottom,
%   bottommost-left, topmost-left.
%   reference: Haralick and Shapiro, Computer and Robot Vision
%   vol I, Addison-Wesley 1992, pp. 62-64.

if ~statsAlreadyComputed.Extrema
    statsAlreadyComputed.Extrema = 1;

    if ndims(L) > 2
        NoNDSupport('Extrema');
        return
    end

    [stats, statsAlreadyComputed] = ...
        ComputePixelList(L,stats,statsAlreadyComputed);

    for k = 1:length(stats)
        pixelList = stats(k).PixelList;
        if (isempty(pixelList))
            stats(k).Extrema = zeros(8,2) + 0.5;
        else
            r = pixelList(:,2);
            c = pixelList(:,1);

            minR = min(r);
            maxR = max(r);
            minC = min(c);
            maxC = max(c);

            minRSet = r == minR;
            maxRSet = r == maxR;
            minCSet = c == minC;
            maxCSet = c == maxC;

            % Points 1 and 2 are on the top row.
            r1 = minR;
            r2 = minR;
            % Find the minimum and maximum column coordinates for
            % top-row pixels.
            tmp = c(minRSet);
            c1 = min(tmp);
            c2 = max(tmp);

            % Points 3 and 4 are on the right column.
            % Find the minimum and maximum row coordinates for
            % right-column pixels.
            tmp = r(maxCSet);
            r3 = min(tmp);
            r4 = max(tmp);
            c3 = maxC;
            c4 = maxC;

            % Points 5 and 6 are on the bottom row.
            r5 = maxR;
            r6 = maxR;
            % Find the minimum and maximum column coordinates for
            % bottom-row pixels.
            tmp = c(maxRSet);
            c5 = max(tmp);
            c6 = min(tmp);

            % Points 7 and 8 are on the left column.
            % Find the minimum and maximum row coordinates for
            % left-column pixels.
            tmp = r(minCSet);
            r7 = max(tmp);
            r8 = min(tmp);
            c7 = minC;
            c8 = minC;

            stats(k).Extrema = [c1-0.5 r1-0.5
                c2+0.5 r2-0.5
                c3+0.5 r3-0.5
                c4+0.5 r4+0.5
                c5+0.5 r5+0.5
                c6-0.5 r6+0.5
                c7-0.5 r7+0.5
                c8-0.5 r8-0.5];
        end
    end

end

%%%
%%% ComputeBoundingBox
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeBoundingBox(L,stats,statsAlreadyComputed)
%   [minC minR width height]; minC and minR end in .5.

if ~statsAlreadyComputed.BoundingBox
    statsAlreadyComputed.BoundingBox = 1;

    [stats, statsAlreadyComputed] = ...
        ComputePixelList(L,stats,statsAlreadyComputed);

    num_dims = ndims(L);

    for k = 1:length(stats)
        list = stats(k).PixelList;
        if (isempty(list))
            stats(k).BoundingBox = [0.5*ones(1,num_dims) zeros(1,num_dims)];
        else
            min_corner = min(list,[],1) - 0.5;
            max_corner = max(list,[],1) + 0.5;
            stats(k).BoundingBox = [min_corner (max_corner - min_corner)];
        end
    end
end

%%%
%%% ComputeSubarrayIdx
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeSubarrayIdx(L,stats,statsAlreadyComputed)
%   Find a cell-array containing indices so that L(idx{:}) extracts the
%   elements of L inside the bounding box.

if ~statsAlreadyComputed.SubarrayIdx
    statsAlreadyComputed.SubarrayIdx = 1;

    [stats, statsAlreadyComputed] = ...
        ComputeBoundingBox(L,stats,statsAlreadyComputed);
    num_dims = ndims(L);
    idx = cell(1,num_dims);
    for k = 1:length(stats)
        boundingBox = stats(k).BoundingBox;
        left = boundingBox(1:(end/2));
        right = boundingBox((1+end/2):end);
        left = left(1,[2 1 3:end]);
        right = right(1,[2 1 3:end]);
        for p = 1:num_dims
            first = left(p) + 0.5;
            last = first + right(p) - 1;
            idx{p} = first:last;
        end
        stats(k).SubarrayIdx = idx;
    end
end

%%%
%%% ComputeEllipseParams
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeEllipseParams(L,stats,statsAlreadyComputed)
%   Find the ellipse that has the same normalized second central moments as the
%   region.  Compute the axes lengths, orientation, and eccentricity of the
%   ellipse.  Ref: Haralick and Shapiro, Computer and Robot Vision vol I,
%   Addison-Wesley 1992, Appendix A.


if ~(statsAlreadyComputed.MajorAxisLength && ...
        statsAlreadyComputed.MinorAxisLength && ...
        statsAlreadyComputed.Orientation && ...
        statsAlreadyComputed.Eccentricity)
    statsAlreadyComputed.MajorAxisLength = 1;
    statsAlreadyComputed.MinorAxisLength = 1;
    statsAlreadyComputed.Eccentricity = 1;
    statsAlreadyComputed.Orientation = 1;

    if ndims(L) > 2
        NoNDSupport({'MajorAxisLength', 'MinorAxisLength', ...
            'Eccentricity', 'Orientation'});
        return
    end

    [stats, statsAlreadyComputed] = ...
        ComputePixelList(L,stats,statsAlreadyComputed);
    [stats, statsAlreadyComputed] = ...
        ComputeCentroid(L,stats,statsAlreadyComputed);

    for k = 1:length(stats)
        list = stats(k).PixelList;
        if (isempty(list))
            stats(k).MajorAxisLength = 0;
            stats(k).MinorAxisLength = 0;
            stats(k).Eccentricity = 0;
            stats(k).Orientation = 0;

        else
            % Assign X and Y variables so that we're measuring orientation
            % counterclockwise from the horizontal axis.

            xbar = stats(k).Centroid(1);
            ybar = stats(k).Centroid(2);

            x = list(:,1) - xbar;
            y = -(list(:,2) - ybar); % This is negative for the
            % orientation calculation (measured in the
            % counter-clockwise direction).

            N = length(x);

            % Calculate normalized second central moments for the region. 1/12 is
            % the normalized second central moment of a pixel with unit length.
            uxx = sum(x.^2)/N + 1/12;
            uyy = sum(y.^2)/N + 1/12;
            uxy = sum(x.*y)/N;

            % Calculate major axis length, minor axis length, and eccentricity.
            common = sqrt((uxx - uyy)^2 + 4*uxy^2);
            stats(k).MajorAxisLength = 2*sqrt(2)*sqrt(uxx + uyy + common);
            stats(k).MinorAxisLength = 2*sqrt(2)*sqrt(uxx + uyy - common);
            stats(k).Eccentricity = 2*sqrt((stats(k).MajorAxisLength/2)^2 - ...
                (stats(k).MinorAxisLength/2)^2) / ...
                stats(k).MajorAxisLength;

            % Calculate orientation.
            if (uyy > uxx)
                num = uyy - uxx + sqrt((uyy - uxx)^2 + 4*uxy^2);
                den = 2*uxy;
            else
                num = 2*uxy;
                den = uxx - uyy + sqrt((uxx - uyy)^2 + 4*uxy^2);
            end
            if (num == 0) && (den == 0)
                stats(k).Orientation = 0;
            else
                stats(k).Orientation = (180/pi) * atan(num/den);
            end
        end
    end

end

%%%
%%% ComputeSolidity
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeSolidity(L,stats,statsAlreadyComputed)
%   Area / ConvexArea

if ~statsAlreadyComputed.Solidity
    statsAlreadyComputed.Solidity = 1;

    if ndims(L) > 2
        NoNDSupport('Solidity');
        return
    end

    [stats, statsAlreadyComputed] = ...
        ComputeArea(stats,statsAlreadyComputed);
    [stats, statsAlreadyComputed] = ...
        ComputeConvexArea(L,stats,statsAlreadyComputed);

    for k = 1:length(stats)
        if (stats(k).ConvexArea == 0)
            stats(k).Solidity = NaN;
        else
            stats(k).Solidity = stats(k).Area / stats(k).ConvexArea;
        end
    end
end

%%%
%%% ComputeExtent
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeExtent(L,stats,statsAlreadyComputed)
%   Area / (BoundingBox(3) * BoundingBox(4))

if ~statsAlreadyComputed.Extent
    statsAlreadyComputed.Extent = 1;

    if ndims(L) > 2
        NoNDSupport('Extent');
        return
    end

    [stats, statsAlreadyComputed] = ...
        ComputeArea(stats,statsAlreadyComputed);
    [stats, statsAlreadyComputed] = ...
        ComputeBoundingBox(L,stats,statsAlreadyComputed);

    for k = 1:length(stats)
        if (stats(k).Area == 0)
            stats(k).Extent = NaN;
        else
            stats(k).Extent = stats(k).Area / prod(stats(k).BoundingBox(3:4));
        end
    end
end

%%%
%%% ComputeImage
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeImage(L,stats,statsAlreadyComputed)
%   Binary image containing "on" pixels corresponding to pixels
%   belonging to the region.  The size of the image corresponds
%   to the size of the bounding box for each region.

if ~statsAlreadyComputed.Image
    statsAlreadyComputed.Image = 1;

    [stats, statsAlreadyComputed] = ...
        ComputeSubarrayIdx(L,stats,statsAlreadyComputed);

    for k = 1:length(stats)
        subarray = L(stats(k).SubarrayIdx{:});
        if ~isempty(subarray)
            stats(k).Image = (subarray == k);
        else
            stats(k).Image = logical(subarray);
        end
    end
end


%%%
%%% ComputePixelList
%%%
function [stats, statsAlreadyComputed] = ...
    ComputePixelList(L,stats,statsAlreadyComputed)
%   A P-by-2 matrix, where P is the number of pixels belonging to
%   the region.  Each row contains the row and column
%   coordinates of a pixel.

if ~statsAlreadyComputed.PixelList
    statsAlreadyComputed.PixelList = 1;

    % Convert the linear indices to subscripts and store
    % the results in the pixel list.  Reverse the order of the first
    % two subscripts to form x-y order.
    In = cell(1,ndims(L));
    for k = 1:length(stats)
        if ~isempty(stats(k).PixelIdxList)
            [In{:}] = ind2sub(size(L), stats(k).PixelIdxList);
            stats(k).PixelList = [In{:}];
            stats(k).PixelList = stats(k).PixelList(:,[2 1 3:end]);
        else
            stats(k).PixelList = zeros(0,ndims(L));
        end
    end
end

%%%
%%% ComputePerimeterCornerPixelList
%%%
function [stats, statsAlreadyComputed] = ...
    ComputePerimeterCornerPixelList(L,stats,statsAlreadyComputed)
%   Find the pixels on the perimeter of the region; make a list
%   of the coordinates of their corners; sort and remove
%   duplicates.

if ~statsAlreadyComputed.PerimeterCornerPixelList
    statsAlreadyComputed.PerimeterCornerPixelList = 1;

    if ndims(L) > 2
        NoNDSupport('PerimeterCornerPixelList');
        return
    end

    [stats, statsAlreadyComputed] = ...
        ComputeImage(L,stats,statsAlreadyComputed);
    [stats, statsAlreadyComputed] = ...
        ComputeBoundingBox(L,stats,statsAlreadyComputed);

    for k = 1:length(stats)
        perimImage = bwmorph(stats(k).Image, 'perim8');
        firstRow = stats(k).BoundingBox(2) + 0.5;
        firstCol = stats(k).BoundingBox(1) + 0.5;
        [r,c] = find(perimImage);
        % Force rectangular empties.
        r = r(:) + firstRow - 1;
        c = c(:) + firstCol - 1;
        rr = [r-.5 ; r    ; r+.5 ; r   ];
        cc = [c    ; c+.5 ; c    ; c-.5];
        stats(k).PerimeterCornerPixelList = [cc rr];
    end

end

%%%
%%% ComputeConvexHull
%%%
function [stats, statsAlreadyComputed] = ...
    ComputeConvexHull(L,stats,statsAlreadyComputed)
%   A P-by-2 array representing the convex hull of the region.
%   The first column contains row coordinates; the second column
%   contains column coordinates.  The resulting polygon goes
%   through pixel corners, not pixel centers.

if ~statsAlreadyComputed.ConvexHull
    statsAlreadyComputed.ConvexHull = 1;

    if ndims(L) > 2
        NoNDSupport('ConvexHull');
        return
    end

    [stats, statsAlreadyComputed] = ...
        ComputePerimeterCornerPixelList(L,stats,statsAlreadyComputed);
    [stats, statsAlreadyComputed] = ...
        ComputeBoundingBox(L,stats,statsAlreadyComputed);

    for k = 1:length(stats)
        list = stats(k).PerimeterCornerPixelList;
        if (isempty(list))
            stats(k).ConvexHull = zeros(0,2);
        else
            rr = list(:,2);
            cc = list(:,1);
            hullIdx = convhull(rr, cc);
            stats(k).ConvexHull = list(hullIdx,:);
        end
    end
end

%%%
%%% ComputePerimeter
%%%
function [stats, statsAlreadyComputed] = ...
    ComputePerimeter(L, stats, statsAlreadyComputed)

if ~statsAlreadyComputed.Perimeter
    statsAlreadyComputed.Perimeter = 1;

    if ndims(L) > 2
        NoNDSupport('ComputePerimeter');
        return
    end

    B = regionboundariesmex(double(L),8);

    for i = 1:length(B)
        boundary = B{i};
        delta = diff(boundary).^2;
        stats(i).Perimeter = sum(sqrt(sum(delta,2)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stats, statsAlreadyComputed] = ...
                ComputePixelValues(I,stats,statsAlreadyComputed)
            
if ~statsAlreadyComputed.PixelValues
    statsAlreadyComputed.PixelValues = 1;
    
    for k = 1:length(stats)
        stats(k).PixelValues = I(stats(k).PixelIdxList);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stats, statsAlreadyComputed] = ...
    ComputeWeightedCentroid(L,I,stats,statsAlreadyComputed)

if ~statsAlreadyComputed.WeightedCentroid
    statsAlreadyComputed.WeightedCentroid = 1;
    
    [stats, statsAlreadyComputed] = ...
        ComputePixelList(L,stats,statsAlreadyComputed);
    
    for k = 1:length(stats)
        sumIntensity = sum(I(stats(k).PixelIdxList));
        numDims = size(stats(k).PixelList,2);
        wc = zeros(1,numDims);
        for n = 1 : numDims
            M = sum(stats(k).PixelList(:,n) .* ...
                double(I(stats(k).PixelIdxList)));
            wc(n) = M / sumIntensity;
        end
        stats(k).WeightedCentroid = wc;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stats, statsAlreadyComputed] = ...
    ComputeMeanIntensity(I,stats,statsAlreadyComputed)

if ~statsAlreadyComputed.MeanIntensity
    statsAlreadyComputed.MeanIntensity = 1;
    
    [stats, statsAlreadyComputed] = ...
        ComputePixelValues(I,stats,statsAlreadyComputed);
    
    for k = 1:length(stats)
        stats(k).MeanIntensity = mean(stats(k).PixelValues);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stats, statsAlreadyComputed] = ...
    ComputeMinIntensity(I,stats,statsAlreadyComputed)

if ~statsAlreadyComputed.MinIntensity
    statsAlreadyComputed.MinIntensity = 1;
    
    [stats, statsAlreadyComputed] = ...
        ComputePixelValues(I,stats,statsAlreadyComputed);
    
    for k = 1:length(stats)
        stats(k).MinIntensity = min(stats(k).PixelValues);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stats, statsAlreadyComputed] = ...
    ComputeMaxIntensity(I,stats,statsAlreadyComputed)

if ~statsAlreadyComputed.MaxIntensity
    statsAlreadyComputed.MaxIntensity = 1;
    
    [stats, statsAlreadyComputed] = ...
        ComputePixelValues(I,stats,statsAlreadyComputed);
    
    for k = 1:length(stats)
        stats(k).MaxIntensity = max(stats(k).PixelValues);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L I reqStats officialStats] = ParseInputs(varargin)

shapeStats = {
    'Area'
    'Centroid'
    'BoundingBox'
    'SubarrayIdx'
    'MajorAxisLength'
    'MinorAxisLength'
    'Eccentricity'
    'Orientation'
    'ConvexHull'
    'ConvexImage'
    'ConvexArea'
    'Image'
    'FilledImage'
    'FilledArea'
    'EulerNumber'
    'Extrema'
    'EquivDiameter'
    'Solidity'
    'Extent'
    'PixelIdxList'
    'PixelList'
    'Perimeter'};

pixelValueStats = {
    'PixelValues'
    'WeightedCentroid'
    'MeanIntensity'
    'MinIntensity'
    'MaxIntensity'};

iptchecknargin(1, inf, nargin, mfilename);

L = varargin{1};
if islogical(L)
    eid = 'Images:regionprops:binaryInput';
    error(eid, ...
        'Use bwlabel(BW) or double(BW) convert binary image to \n%s', ...
        'a label matrix before calling regionprops.');
end
iptcheckinput(L, {'numeric'}, ...
    {'real', 'integer', 'nonnegative','nonsparse'}, ...
    mfilename, 'L', 1);

basicStats = {'Area'
    'Centroid'
    'BoundingBox'};

I = [];
officialStats = shapeStats;
if nargin == 1
    %REGIONPROPS(L)
    reqStats = basicStats;
    return;
elseif isnumeric(varargin{2}) || islogical(varargin{2})
    %REGIONPROPS(L,I) or REGIONPROPS(L,I,PROPERTIES)
    I = varargin{2};
    iptcheckinput(I, {'numeric','logical'},{}, mfilename, 'I', 2);
    assert(isequal(size(L),size(I)), ...
        'Images:regionprops:sizeMismatch', ...
        'I must have the same size as L.');
    
    officialStats = [shapeStats;pixelValueStats];
    if nargin == 2
        %REGIONPROPS(L,I)
        reqStats = basicStats;
        return;
    else
        %REGIONPROPS(L,I,PROPERTIES)
        startIdxForProp = 3;
        reqStats = getPropsFromInput(startIdxForProp, ...
            officialStats, pixelValueStats, basicStats, varargin{:});
    end
else
    %REGIONPROPS(L,PROPERTIES)
    startIdxForProp = 2;
    reqStats = getPropsFromInput(startIdxForProp, ...
        officialStats, pixelValueStats, basicStats, varargin{:});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [reqStats,officialStats] = getPropsFromInput(startIdx, officialStats, ...
    pixelValueStats, basicStats, varargin)

if iscell(varargin{startIdx})
    %REGIONPROPS(...,PROPERTIES)
    propList = varargin{startIdx};
elseif strcmpi(varargin{startIdx}, 'all')
    %REGIONPROPS(...,'all')
    reqStats = officialStats;
    return;
elseif strcmpi(varargin{startIdx}, 'basic')
    %REGIONPROPS(...,'basic')
    reqStats = basicStats;
    return;
else
    %REGIONPROPS(...,PROP1,PROP2,..)
    propList = varargin(startIdx:end);
end

numProps = length(propList);
reqStats = cell(1, numProps);
for k = 1 : numProps
    if ischar(propList{k})
        noGrayscaleImageAsInput = startIdx == 2;
        if noGrayscaleImageAsInput
            % This code block exists so that regionprops can throw a more
            % meaningful error message if the user want a pixel value based 
            % measurement but only specifies a label matrix as an input.
            tempStats = [officialStats; pixelValueStats];
            prop = iptcheckstrs(propList{k}, tempStats, mfilename, ...
                'PROPERTIES', k);
            if any(strcmp(prop,pixelValueStats))
                eid = sprintf('Images:%s:needsGrayscaleImage', mfilename);
                error(eid, 'REGIONPROPS needs I as an %s ''%s''.', ...
                    'input to calculate', prop);
            end
        else
            prop = iptcheckstrs(propList{k}, officialStats, mfilename, ...
                'PROPERTIES', k);
        end
        reqStats{k} = prop;
    else
        eid = sprintf('Images:%s:invalidType', mfilename);
        error(eid, 'PROPERTY must be a string.');
    end
end


%%%
%%% NoNDSupport
%%%
function NoNDSupport(str)
%   Issue a warning message about lack of N-D support for a given
%   measurement or measurements.

wid = sprintf('Images:%s:measurementNotForN-D',mfilename);

if iscell(str)
    warn_str = sprintf('%s: %s ', ...
        'These measurements are not supported if ndims(L) > 2.', ...
        sprintf('%s ', str{:}));
else
    warn_str = sprintf('%s: %s', ...
        'This measurement is not supported if ndims(L) > 2.', ...
        str);
end

warning(wid,'%s',warn_str);

%%%
%%% PreprocessRequestedStats
%%%
function requestedStats = PreprocessRequestedStats(requestedStats)
%   Remove any requested stats that are not supported for N-D input
%   and issue an appropriate warning.

no_nd_measurements = {'MajorAxisLength'
    'MinorAxisLength'
    'Eccentricity'
    'Orientation'
    'ConvexHull'
    'ConvexImage'
    'ConvexArea'
    'EulerNumber'
    'Extrema'
    'EquivDiameter'
    'Solidity'
    'Extent'
    'Perimeter'};

bad_stats = find(ismember(requestedStats, no_nd_measurements));
if ~isempty(bad_stats)
    NoNDSupport(requestedStats(bad_stats));
end

requestedStats(bad_stats) = [];
