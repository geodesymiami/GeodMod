function []=make_kmz(data,rgbim,seed);
% make_kmz     - creates a kmz file containing a jpeg image and kml file 
%
% usage: []=make_kmz(igram,rgbim);
%
%
%  Scott Baker, May 2007
%  
%   Modified: N. Gourmelen - Feb 2008
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dir_out


%[junk,junk,fname]=GenerateNamesFromDataStructure(data);
filename = ['timeseries'] ;
if ~isfield(data,'date_years')  data.date_years = 1 ;  end
if exist('overlays') ~= 7;  unix('mkdir overlays');  end


%% Check min and max bounds for all dataset

all_max = -1e100 ;  all_min = 1e100 ;

for ni=1:length(data)    
    [mmax,mmin] = maxmin(data(ni).data);
    if (mmax > all_max);  all_max = mmax ;  end
    if (mmin < all_min);  all_min = mmin ; end
end

bound = max([abs(all_min) abs(all_max)]);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Get the Coordinates %%%%%%%%%
% label vectors for 'Coord'
xl          = linspace(data(1).x_first,data(1).x_first+data(1).x_step*size(data(1).data,2),size(data(1).data,2)) ;
yl          = linspace(data(1).y_first,data(1).y_first+data(1).y_step*size(data(1).data,1),size(data(1).data,1)) ;
bounds      = [xl(1) data(1).x_step xl(size(xl,2)) yl(1) data(1).y_step yl(size(yl,2))]                          ;
north_coord = bounds(4) ;
south_coord = bounds(6) ;
east_coord  = bounds(3) ;
west_coord  = bounds(1) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Create the KML string %%%%%%%%%%
header1 = ['<?xml version="1.0" encoding="UTF-8"?>\n',...
    '<kml xmlns="http://earth.google.com/kml/2.1">\n',...
    '<Document>\n',...
    '<name>\n',...
    filename,...
    '\n</name>\n',...
    '<description>\n',...
    'Image overaly created with geodmod',...
    '</description>']

init_string = [];

for ni=1:length(data)
    
    [yr,mm,dd] = jd2cal(yr2jd(data(ni).date_years));

    if ( ni < length(data) )  [yr_end,mm_end,dd_end] = jd2cal(yr2jd(data(ni+1).date_years)) ;
    else                      [yr_end,mm_end,dd_end] = jd2cal(yr2jd(data(ni).date_years))   ;
    end
    
    date_start = [sprintf('%4.0f',yr),'-',sprintf('%02.0f',mm),'-',sprintf('%02.0f',dd)]             ;
    date_end   = [sprintf('%4.0f',yr_end),'-',sprintf('%02.0f',mm_end),'-',sprintf('%02.0f',dd_end)] ;
    
    image_file = ['overlays/',date_start,'.png'];  % changed image filename to be the date of image, sbaker:Feb2008
    colfact=1;  % Factor to shift the color scale
    %[max,min] = maxmin(data(ni).data)
    tmp_data = data(ni).data * colfact ;  %[mmax,mmin] = maxmin(tmp_data) ;  
    tmp_data = (tmp_data-mmin)/(2*bound); 
    mean_data = mean(mean(tmp_data(isfinite(tmp_data)))) ; 

    if ~isfield(data,'seed') seeder = 0 ;  else  seeder = tmp_data(data(ni).seed(1),data(ni).seed(2)) ;  end

    tmp_data  = round((tmp_data+(0.5-seeder))*64) ; 
    tmp_data(tmp_data>64)=64;
    [nan_ind] = find(isnan(tmp_data)==1) ;  %tmp_data(nan_ind) = 0 ;
    rgbim_im = ind2rgb(tmp_data,rgbim);  mask = ones(size(tmp_data)) ;  mask(nan_ind) = 0 ;
    
    save png_maker_files
    imwrite(rgbim_im,image_file,'PNG','Alpha',mask,'BitDepth',16)

    tmp = ['<GroundOverlay>\n',...
        '<name>\n',...
        image_file,...
        '</name>\n',...
        '<Icon><href>\n',...
        image_file,...
        '\n</href></Icon>\n',...
        '<LatLonBox>\n',...
        sprintf('\t<north>%.6f</north>',north_coord),...
        sprintf('\n\t<south>%.6f</south>', south_coord),...
        sprintf('\n\t<east>%.6f</east>', east_coord),...
        sprintf('\n\t<west>%.6f</west>', west_coord),...
        '\n</LatLonBox>\n',...        
        '<TimeSpan><begin>',date_start,'</begin><end>',date_end,'</end></TimeSpan>',...
        '</GroundOverlay>\n'];

    init_string = strcat(init_string,tmp);

end

footer1 = ['</Document>\n','</kml>\n'] ;

kml_string = strcat(header1,init_string,footer1) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Create the KMZ File %%%%%%%%%%
%imwrite(rgbim,'image_overlay.jpg','JPEG')        % create the image overlay
kml_file = 'doc.kml';                             % create the kml file
%image_file = 'image_overlay.jpg';
filename = strcat(filename,'.kmz');
% Write the contents of kml_string to kml_file
fid = fopen( kml_file, 'wt');
fprintf(fid, kml_string );
fclose(fid);
% Create and zip the KMZ File


zip(filename,{kml_file,'overlays'});                  % zip the kml file and image jpeg
zipped_file = strcat(filename,'.zip');            % the zip command appends .zip to the file name which we don't want
movefile(zipped_file,filename);                   % move the zipped file to the name ending with .kmz
%delete(kml_file,image_file);                      % delete the kml and image file to save space since they are in the .kmz file now

logmessage (sprintf('Created GoogleEarth KMZ File: %s',filename))




