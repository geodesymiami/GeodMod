function [data_type]=Data2Igram(opt)
%Data2Igram           - reads roi_pac-processed interferograms and saves igram structure
%
%  usage:  Data2Igram(opt)
%
%          Options given as opt.dir_in, etc.
%
%          'dir_in'            Input directory (contains e.g. RsatA3)       [default 'off']
%          'pattern'           filematch                                    [default 'PO*/geo*[0-9].unw']
%          'igram_name'        filename for igram file                      [default is obtained with extract_name_fromSOdir]
%          'rating'            'on' or 'off' for use of igram_rating.log    [default 'off']
%          'unw_thresh'        unwrap threshold for LoadData                [default 0.8]
%          'incangleFile'file containing inc angle and azinuth        [default 'geo_incidence.unw']
%          'corfile'           coherence to mask snaphu-unwrapped data      [default 'off']
%                              if 'off' the regular geo*cor is used
%          'cormask_thresh'    threshold for coherence mask                 [default 0.2]
%          'cormask_area_open' minimum connectd pixel in coherence mask     [default 10]
%          'Unit'              unit of data                                 [default 'radian']
%          'plotdataopt'       plot options                                 [default 'off' ]
%          'ambiguityOffset'   phase ambiguity area extracted from corr                 [default 'off' ]
%          'SubtractModelFile'                                              [default 'off' ]
%                        
%  Saves data as igram structure.
%  If corfile is given it is used 
%  as a mask. Program could be easily 
%  changed to accept another mask.
%  Removes mean after applying mask.  
%  
%  Part of the TimeSeries suite
%  V1.0  Falk Amelung, September 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ls -1 */geo*[0-9].cor > ! cor.list  ; meanimage cor.list 2551 4786 mean.cor ; cp `head -1 cor.list`.rsc mean.cor.rsc ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultopt=struct(                                                         ...
        'DoIt'               ,        'on'                   ,            ...
        'in_name'            ,        'off'                   ,            ...
        'out_name'           ,        'off'                   ,            ...
        'pattern'            ,        ['PO*' filesep 'geo*[0-9].unw']     ,            ...
        'igram_name'         ,        'off'                   ,            ...
        'rating'             ,        'off'                   ,            ...
        'unw_thresh'         ,         0.99                   ,            ...
        'corfile'            ,        'off'                   ,            ...
        'cormask_thresh'     ,         0.2                    ,            ...
        'cormask_area_open'  ,         10                     ,            ...
        'Unit'               ,        'radian'                ,            ...
        'SubtractModelFile'  ,        'off'                   ,            ...
        'ambiguityOffset'    ,        'off'                   ,            ...
        'ambiguityOffsetFile',        'off'                   ,            ...
        'plotdataopt'        ,        'off'      )            ;
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],defaultopt); end ;
[opt]=process_defaultoptions(opt,defaultopt);  %display(opt)
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end
if  ~DoIt                    return; end
if  CheckInOut('',out_name)  tmp = load(out_name) ;  data_type = tmp.(char(fieldnames(tmp))).Unit ;  return; end

dir_in=in_name ;
if exist(dir_in,'dir') filematch=fullfile(dir_in,pattern); else filematch=dir_in; end           % 10/2009: remove once Pattern option is eliminated

if ~igram_name
 [junk,junk,satbeam]=extract_ProjectName(dir_in);
 igram_name         = [ 'igram_' satbeam  '.mat'] ;
end

%
% load interferograms
%
   loadopt.Unit=Unit; 

   [igram,datelist,N_timevector]=LoadData(filematch,loadopt);
   
%
% load incidence angle file if it exists (geo_incidence.unw) (FA 1/2011)
%

    [pathstr, name, ext] = fileparts(filematch);
    filematch_incangle   = [pathstr filesep 'geo_incidence.unw'];
      
    if (exist(filematch_incangle,'file')==2 || ~isempty(findstr(filematch_incangle,'*')))    % Load only if file 
       loadopt.loadAmplitude= 1;
       loadopt.dataMean     = 1;
       [tmp]                = LoadData(filematch_incangle,loadopt);
       %igram(1).inc1        = tmp(1).amp;
       igram(1).azimuth     = tmp(1).dataMean;
       igram(1).incangle    = mean(reshape([tmp(:).amp],size(tmp(1).amp,1),size(tmp(1).amp,2),[]),3); %  calculates mean incangle arrays  
    end
   %
   % load coherence and apply threshold
   %
      if corfile
         corfile=fullfile(dir_in,corfile);
         coropt.Unit='unity';
         coropt.DataType='Coherence';
         [correlation]=LoadData(corfile,loadopt);
      else
         [pathstr,name,ext]   = fileparts(filematch);
         corfilematch=fullfile(pathstr,[name '.cor']);
         [correlation,datelist,N_timevector]=LoadData(corfilematch,loadopt);
      end

   %
   % theshold coherence and remove small objects from mask
   %

         for i=1:length(correlation)

             if exist('cormask_thresh_value')                          % FA 1/11 full switch to in_list
                thresh = cormask_thresh_value;
             elseif isstruct(cormask_thresh)                         % FA 1/11 for backward compatibility (need to remove)
                if ~isfield(cormask_thresh,satbeam) errordlg(sprintf('cormask.%s not specified, -- exiting',satbeam)); end
                thresh = cormask_thresh.(satbeam);
             else 
                thresh = cormask_thresh;
             end

             correlation(i).data(find(correlation(i).data< thresh)) = nan;    % thresholding    
             correlation(i).data(find(correlation(i).data>=thresh)) = 1;

             if cormask_area_open                                     % remove small objects
                cormask = correlation(i).data;
                cormask(find(correlation(i).data==1))    =0;           % assign 0,1 to high,low coherence pixel so  
                cormask(find(isnan(correlation(i).data)))=1;           % that bwareaopen works
                tmp=bwareaopen(cormask,cormask_area_open)*1.0 ;
                cormask(find(tmp==1))  =NaN;                           % assign NaN to low  coherence pixel
                cormask(find(tmp==0))  =1;                             %         1  to high coherence pixel
             end
             
             correlation(i).data = cormask ;
         end

   %
   % apply coherence mask
   %
         for i=1:length(igram)
             if corfile j=1 ; else j=i ; end
             
             %%% Tini 6/2010: we may want to give all the ambiguity stuff
             %%% in a structure so that multiple ambiguites can be given
             [igram(i)]=add_mask_to_igram(igram(i),correlation(j).data,loadopt,opt);
         end

   %
   % add rating
   %
   if rating
       if islogical(rating)
          ratingfile=fullfile(dir_in,'igram_rating.log');
       else
          ratingfile=fullfile(dir_in,rating);
       end
       if exist(ratingfile,'file')   
           [igram]=add_rating_to_igram(igram,ratingfile);  
        else
            error('user error: ratingfile %s does not exist -- exiting',ratingfile)
       end
   end

   %[igram]=remove_mean(igram); %FA 12/18: Commented out as it looks unnecessary

%
% add DataSet to igram structure
%
   [junk,junk,DataSet]=extract_ProjectName(dir_in);
   for i=1:length(igram)
       igram(i).DataSet  = DataSet;
   end

   %plotdataopt.CLim='Centered';
   if isfield(plotdataopt,'CLim') plotdataopt=rmfield(plotdataopt,'CLim') ; end
   logplot('PlotData_wrapper',out_name,igram,plotdataopt);

   data_type = igram(1).Unit;
 
   
   
%--------------------------------------------------------------------------
% Calculate and Subtract a forward model from each Igram (e.g. Interseismic signal): 
% Tini 09/2010
 if (opt.SubtractModelFile)   
      for i=1:length(igram)
      [igram(i)]=SubtractForwardModel(igram(i),opt);
      end
end
%--------------------------------------------------------------------------


% save data
%
  save(out_name,'igram')
