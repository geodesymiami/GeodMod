function [indata] = Geo2Radar(indata,opt);

%
% Add radar coordinates to the dtaa structure
%
%            [outdata]=Read_IREA(indata,opt)
%
%
%       indata : Input file
%
%
% opt contains 
%
%   'gps'           : Location of GPS files
%
%   'gpsradcoor'    : name of the file containing the date list
%
%
%   N. Gourmelen, April 2009
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultopt=struct(                                  ...
                    'gps'                            ,             'off'           ,           ...
                    'radcoor'                        ,             'off'           )            ;


if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],[]); end ;
[opt]=process_defaultoptions(opt,defaultopt);  display(opt)
f=fieldnames(opt) ;
for i=1:length(f)
    eval([char(f{i}) '= opt.(f{i}) ;' ]) ;
end

if gpsradcoor
    [filepath,file,ext] = fileparts(gps)        ;
    fid = fopen([filepath,'/',gpsradcoor],'r')  ;
    C   = textscan(fid,'%s%f%f')                ;
    
    for ni = 1:length(indata)
        tmpi = find(strcmp(C{1},indata(ni).GPS_station)==1) ;
        jC = C{2} ;  iC = C{3} ;
        if ~isempty(tmpi)
            indata(ni).nj = jC(tmpi) ;  indata(ni).ni = iC(tmpi) ;
        else
            indata(ni).nj = -1 ;  indata(ni).ni = -1 ;
        end
    end
end
    
    
