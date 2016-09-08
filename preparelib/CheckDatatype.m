function [data_type]=CheckDatatype(opt)

%   Returns Unit from a dry run (could be faster for IREA inputs)
%
%  usage:  [data_type]=CheckDatatype(opt)
%
%          Options given as opt.dir_in, etc.
%
%          'dir_in'            Input directory (contains e.g. RsatA3)       [default 'off']
%          'pattern'           filematch                                    [default 'PO*/geo*[0-9].unw']
%          'igram_name'        filename for igram file
%          [default is obtained with extract_name_fromSOdir]
%          [default 'off']
%          'Unit'              unit of data                                 [default 'radian']
%                        
%  
%  Part of the TimeSeries suite
%  V1.0  Noel Gourmelen - Feb 2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultopt=struct(                                                         ...
        'DoIt'               ,        'on'                   ,            ...
        'in_name'            ,        'off'                   ,            ...
        'out_name'           ,        'off'                   ,            ...
        'pattern'            ,        ['PO*' filesep 'geo*[0-9].unw']     ,            ...
        'igram_name'         ,        'off'                   ,            ...
        'Unit'               ,         'radian'               )             ;
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],defaultopt); end ;
[opt]=process_defaultoptions(opt,defaultopt);  %display(opt)
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end
if  ~DoIt                    return; end

if exist(out_name,'file')
    S         = load(out_name);
    data_type = S.(char(fieldnames(S))).Unit;
    logmessage(sprintf('data_type found in: %s',out_name));
else
    dir_in=in_name ;
    if exist(dir_in,'dir') filematch=fullfile(dir_in,pattern); else filematch=dir_in; end
    %
    % load interferograms
    %
    loadopt.Unit = Unit ;  loadopt.DryRun = 'on' ;  [data_type] = LoadData(filematch,loadopt) ;
end

   
