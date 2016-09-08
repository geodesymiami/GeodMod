function []=Igram2Motion(opt)
%Igram2Motion             - converts igram structure into rate using different stacking options
%
%  usage:  Igram2Motion(opt)
%
%          optional arguments (given as 'dir_in','RAID6/SO/HawaiiRsatA3,'pattern','PO*/geo*[0-9].unw',...)
%
%          'igram_name'         Input directory (contains e.g. RsatA3)       [default 'off'  ]
%          'dir_out'            Output directory                             [default wdir   ]
%          'outfile'            filename for igram file                      [default is obtaining from dir_in if given]
%          'Rating'             Maximum rate included                        [default 9999   ]
%          'tforstack'          time period for averaging rate               [default 'off   ]
%                               ('20020512') or ('20020512':20020531')                            
%                               (second date not yet implemented)                            
%          'RateMethod'         method to generate rate (SimpleAverage,Rand) [default 'Rand' ]
%          'RandState'          Seed for random generator ('off' for 1)      [default  1     ]
%          'NumStacks'          Number of random stacks                      [default  10    ]
%          'NumImRepeat'        How many times the same image                [default  1     ]
%          'NumIntMin'          Minimum number of Interferograms in stack    [default  2     ]
%          'TimeSpanMin'        Minimum time span between images(in yr)      [default  0     ]
%          'TimeSpanMax'        Maximum time span between images(in yr)      [default 9999   ]
%          'CSpaceLim'          Window for calculating colorscale [x y]      [default all    ]          
%          'CLim'               Colorscale bounds                            [default min/max]
%          'StackSortMethod'Method to sort the stack                 [default 'MeanRating']
%                           (MeanRating,TotalTime)
%          'MedFiltSize'        final median filtersize after stacking       [default min/max]
%          'Unit'               unit of data, mm/yr,cm/yr,m/yr               [default 'm/yr' ]
%          'plotdataopt'       plot options                                 [default 'off' ]
% 
%  Saves data as igram structure.
%  If corfile is given it is used 
%  as a mask. Program could be easily 
%  changed to accept another mask.
%  Removes mean after applying mask.  
%  
%  Part of the TimeSeries suite
%  V1.0  Noel Gourmelen, Falk Amelung, September 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  set default options. Process options fom input or from file *.min
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultopt=struct(                                                         ...
        'DoIt'               ,        'on'                   ,            ...
        'in_name'            ,        'off'                   ,            ...
        'out_name'           ,        'off'                   ,            ...
        'dir_out'            ,        'off'                   ,            ...
        'Rating'             ,        9999                    ,            ...
        'RateMethod'         ,        'Rand'                  ,            ...
        'RandState'          ,        0                       ,            ...
        'MedFiltSize'        ,        'off'                   ,            ...
        'NumStacks'          ,        10                      ,            ...
        'NumImRepeat'        ,        9999                    ,            ...
        'NumIntMin'          ,        2                       ,            ...
        'TimeSpanMin'        ,        0                       ,            ...
        'TimeSpanMax'        ,        9999                    ,            ...
        'tforstack'          ,        'off'            ,            ...
        'StackSortMethod'    ,        'MeanRating'            ,            ...
        'CSpaceLim'          ,        'off'                   ,            ...
        'Unit'               ,        'm/yr'                  ,            ...
        'plotdataopt'        ,        'off'      )            ;
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],defaultopt); end ;
[opt] = process_defaultoptions(opt,defaultopt);  %display(opt)
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end
logmessage(sprintf('[]=%s(%s)',mfilename,inputname(1)));
if  ~DoIt                         return; end
if  CheckInOut(in_name,out_name)  return; end
%
% Load data and adjust igram structure and input option 
%
load(in_name,'igram')

%  add tmin and tmax to igram structure if given
if tforstack & ~isempty(tforstack)
    [igram]=add_tforstack_to_igram(igram,tforstack);
end

if length(igram)==1 
   RateMethod='SimpleAverage';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Noel, all the following should be in one function. The function should understand the options RateMethod and Unit
% It also should understand
% It should make 2 plotts: stacks_RsatA4.png and stacks_RsatA4_mul.png The first is the selected stack, The second shows all
% stacks.

if ~strcmp(igram(1).Unit,'radian') error('unit of igram needs to be radian - existing ] \n'); end   

wavelength = extract_hardwired_satparameters(igram(1).DataSet,'wavelength')

if       strcmp(Unit,'m/yr')   Fac=wavelength/4/pi * 1    ;   rate_flag=true;    %factor for conversion into m/yr
  elseif strcmp(Unit,'cm/yr')  Fac=wavelength/4/pi * 100  ;   rate_flag=true;    %factor for conversion into cm/yr
  elseif strcmp(Unit,'mm/yr')  Fac=wavelength/4/pi * 1000 ;   rate_flag=true;    %factor for conversion into mm/yr
  elseif strcmp(Unit,'m')      Fac=wavelength/4/pi * 1    ;   rate_flag=false;   %factor for conversion into m
end

% Noel, we need a function [Fac]=UnitConversion(unit_in,unit_out) that calculates Fac depending on the Units
% it also should understands unit_in=[].

switch RateMethod
case {'SimpleAverage'}
     [S]=SimpleAverage(igram,'Fac',Fac,'PlotStack','off','CLim',[-0.02 0.02],'Rate',rate_flag) ;
case {'Rand'}

       rand('state',RandState);

       oldstruc.NbStck       = NumStacks ;
       oldstruc.rating       = Rating ;
       oldstruc.repeatNb     = NumImRepeat;
       oldstruc.LowTmpThresh = TimeSpanMin;
       oldstruc.HighTmpThresh= TimeSpanMax;
       oldstruc.MinNbInt     = NumIntMin;
       %oldstruc.MultiFact    = Fac;
       %oldstruc.CSpaceLim    = CSpaceLim;

       %[S]=StackRandFalk(igram,'NbStck',NumStacks,'rating',Rating,'DoPlotIgram','off','PlotStack','on','CLim',[-0.02 0.02]) ;

       [S]=StackRandFalk(igram,'NbStck',NumStacks,'rating',Rating,'Fac',Fac,'StackSortMethod',StackSortMethod,'Rate',rate_flag) ;
end

if MedFiltSize 
   logmessage(sprintf('median filtering, MedFiltSize: %d, can be switched off',MedFiltSize))
   %logmessage(sprintf('search more efficient function than nanmedfilt2 (matlab's medfilt2 is fast but does not ignore NaNs'))
   tmp=nanmedfilt2(S(1).data,MedFiltSize,MedFiltSize,1000); S(1).data=tmp;
end

tmpDates              =  char(S(1).IgramDates{:});

motion.data           = S(1).data;
motion.x_first        = S(1).x_first;
motion.y_first        = S(1).y_first;
motion.x_step         = S(1).x_step;
motion.y_step         = S(1).y_step;
motion.x_unit         = S(1).x_unit;
motion.DataSet        = igram(1).DataSet;
motion.IgramIndices   = S(1).IgramIndices ; 
motion.TotalTime      = S(1).TotalTime/length(S(1).IgramIndices);
motion.StartDate      = mean(YYYYmmDD2DecimalYear(tmpDates(:,1:8)));
motion.StackTotalTime = S(1).TotalTime;
motion.IgramDates     = S(1).IgramDates ; 
motion.Unit           = Unit ;

if isfield(igram,'incangle')  motion.incangle  = igram(1).incangle; end      % FA 1/2011: To be correct we should average over incidence angle
if isfield(igram,'azimuth')   motion.azimuth   = igram(1).azimuth;  end      % FA 1/2011: needed to calculate radarlook at each pixel

clear igram opt          % free some memory

%plotdataopt.googleearthopt.DoIt = 'on' ;

logplot('PlotData',out_name,motion,plotdataopt);

% plot the Stacks if S contains them
  if length(S) > 1
     %[pathstr,name,ext]     = fileparts(out_name) ;   
     %plotname=[pathstr filesep 'stacks_' motion.DataSet];           % mul stands for multiple
     
     %logplot('PlotData_wrapper',out_name,igram,plotdataopt);
     % Commented out by Tini because of unknown variable plotname 10/2010
%    logplot('plot_igram',plotname,S,'FigName','Stacks','Coord','off');
     
      [S(:).Unit]    = deal(Unit);
      [S(:).DataSet] = deal(motion.DataSet);
     
      logplot('PlotData_wrapper',out_name,S(:),plotdataopt);
     %logplot('PlotData_wrapper',out_name,S,'FigName','Stacks','Coord','off');

  end

save(out_name,'motion')

