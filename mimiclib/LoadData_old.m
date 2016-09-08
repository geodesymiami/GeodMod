function [igram,datelist,N_timevector]=LoadData(filematch,opt)
%LoadData               - Loads Data into igram structure
%
%   Input options:
%
%   opt.subset    subset to be loaded (e.g.: [250 250 512 512]              [default : entire file]
%   opt.unw_thresh    unwrap threshold, percentage needs to be above this value [default: 0.2 ]
%   opt.remove    
%   opt.Unit                 Unit of data ([radian,mm,cm,m,unit])                [default 'radian']
%                             note that sign will be changed except for 'unit' (use for coherence)
%   opt.cormask_thresh       assigns 0 below and 1 above threshold value
%   opt.cormask_area_open    removes to small areas (useful to remove individual pixel with high coherence)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ATTENTION:  PROGRAM NEEDS TO BE CONVERTED INTO THE TYPICAL OPTION STRUCTURE BUT NO TIME AS OF SEP 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Nov 14 2005, FA, converting some striungs to numbers
%   Sep 2006,    FA, added options opt.cormask_thresh, opt.cormask_area_open for making coherence masks
%   18 Sep 2006, FA, replaced option opt.fac by opt.Unit
%   06 Oct 2006, FA, for Unit='radian' now inverting sign of data so that negative phase (LOS decrease) is positive (uplift)
%                    This corresponds to a change of the coordinate system from satellite-based to ground-based.
%   14 Aug 2007  FA, replaced ReadRscfile by ReadKeywordfile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultopt=struct(                                                         ...
        'subset'             ,        'off'                   ,            ...
        'remove'             ,        'off'                   ,            ...
        'unw_thresh'         ,         0.99                    ,            ...
        'Unit'               ,        'radian'                ,            ...
        'DataType'           ,        'off'                ,            ...
        'cordata'            ,        'off'                   ,            ...
        'cormask_thresh'     ,         0.2                    ,            ...
        'cormask_area_open'  ,         10        )            ;
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],defaultopt); end ;
[opt]=process_defaultoptions(opt,defaultopt);  %display(opt)
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval(sprintf('system(''ls %s > tmp_filelist'');',filematch)) ;
[filelines]=textread('tmp_filelist','%s','delimiter','\n');
logmessage('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% start loop to load data %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(filelines)
   file=filelines{i};

   logmessage(sprintf('Loading file %d %s',i, file))

   rscfile=strcat(file,'.rsc');
   S=ReadKeywordfile(rscfile,' ');
   width   = S.WIDTH  ;
   x_first = S.X_FIRST;
   y_first = S.Y_FIRST;
   x_step  = S.X_STEP ;
   y_step  = S.Y_STEP ;
   x_unit  = S.X_UNIT ;

   d=[];
   if isfield(S,'DATE12')
       date12=S.DATE12 ;
       d=findstr(date12,'-');
   end


   if subset
      [a,p]=readhgt(file,width,opt.subset);
   else
      [a,p]=readhgt(file,width);
   end
   p(isnan(a))=nan ;
   p(find(a==0))=nan;
   delta_a=0.01;
   p(find(a>0-delta_a & a<0+delta_a))=nan;
   
   if subset
      x_first=x_first+opt.subset(1)*x_step;
      y_first=y_first+opt.subset(2)*y_step;
   end

   igram(i).data=p;               
   igram(i).x_first=x_first;
   igram(i).y_first=y_first;
   igram(i).x_step=x_step;
   igram(i).y_step=y_step;
   igram(i).x_unit=x_unit;

   %
   % add date information to igram structure if infile is interferogram (d not empty)
   %
   if ~isempty(d) 
      date1=date12(1:d-1) ; 
      date2=date12(d+1:end);
      sd1=num2str(date1) ;  sd2=num2str(date2);
      mjd1=date2j(str2num(sd1(1:2)),str2num(sd1(3:4)),str2num(sd1(5:6)));
      mjd2=date2j(str2num(sd2(1:2)),str2num(sd2(3:4)),str2num(sd2(5:6)));
   
      % change sign of that first acquisition is first image if necessary
      if mjd1 > mjd2
         p=-p        ;
         qtmp=date1  ; date1=date2 ; date2=qtmp ; 
         qtmp=mjd1   ;  mjd1=mjd2  ; mjd2 =qtmp ;
       end

      if (strcmp(date1(1),'9')==1)   date1=strcat('19',date1) ; end    %change march 05
      if (strcmp(date1(1),'0')==1)   date1=strcat('20',date1) ; end
      if (strcmp(date2(1),'9')==1)   date2=strcat('19',date2) ; end
      if (strcmp(date2(1),'0')==1)   date2=strcat('20',date2) ; end

      igram(i).t1=mjd1;
      igram(i).t2=mjd2;
      igram(i).delt=mjd2-mjd1;
      igram(i).date1=date1;
      igram(i).date2=date2;
   end

   %
   % fill cordata with ones if we are reading a Dem
   %  this should be done in PreparePlot but the matrix size is not available because of the DemAreaEqualsSubset option
   if strcmp(DataType,'Dem') | strcmp(DataType,'Coherence')
      cordata=ones(size(igram(1).data)) ;
      opt.cordata=cordata ;
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%% process coherence data/file   %%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %  read coherence file if no cordata passed via opt.cordata
   %
   if islogical(opt.cordata) && ~opt.cordata   
      [pathstr,name,ext]   = fileparts(file);
      corfile=fullfile(pathstr,[name '.cor']);
      if isfield(opt,'subset')
          [ca,cordata]=readhgt(corfile,width,opt.subset); 
      else
          [ca,cordata]=readhgt(corfile,width,opt.subset); 
      end
   end

   %
   % apply coherence theshold  
   %
   if isfield(opt,'cormask_thresh')
      cormask_thresh=opt.cormask_thresh ;
      cordata(find(cordata< cormask_thresh))=nan;
      cordata(find(cordata>=cormask_thresh))=1;
   end
   %
   % remove small objects 
   %
    if isfield(opt,'cormask_area_open')
       area_open=opt.cormask_area_open ;
       cor_bw=im2bw(cordata) ;
       tmp=bwareaopen(cor_bw,area_open)*1.0 ;
       tmp(find(tmp==0.0))=nan;
       cordata=tmp ;
    end
    [igram(i)]=add_mask_to_igram(igram(i),cordata,opt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logmessage(sprintf('igrams read in: %d',length(igram)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  remove interferograms with too much non-unwrapped data

if unw_thresh
   [igram,N_igrams]=RemoveBadIgrams(igram,unw_thresh);
   logmessage(sprintf('good igrams   : %d',N_igrams))
else
   N_igrams=length(igram);
end

%
%  remove interferograms using remove list
%
if remove
   [igram,N_igrams]=RemoveIgrams(igram,opt.remove);
   logmessage(sprintf('igrams after removing : %d',N_igrams))
end

%
%  convert phase into unit desired
%
if ~isfield(opt,'Unit')
    opt.Unit='radian' ;
end
if         strcmp(opt.Unit,'radian') fac=-1    ;                               %factor for conversion into radian
    elseif strcmp(opt.Unit,'mm')     fac=-28.3/2/pi    ;                       %factor for conversion into mm
    elseif strcmp(opt.Unit,'cm')     fac=-2.83/2/pi    ;                       %factor for conversion into cm
    elseif strcmp(opt.Unit,'m')      fac=-0.0283/2/pi  ;                       %factor for conversion into m
    elseif strcmp(opt.Unit,'unity')   fac=1             ;                       %no change (coherence, height)
end

for i=1:length(igram)
   igram(i).data=igram(i).data*fac;
   igram(i).Unit=opt.Unit;
end

%
%  and sort interferograms
%
if ~isempty(d)                  %sort according to tome only if date was given
   [igram,N_igrams]=SortIgram(igram);
   %  set time of first image to zero and create datelist
   for i=1:N_igrams
       t(i,:)=[igram(i).t1 igram(i).t2];
   end
   firsttime=min(min(t));
   for i=1:N_igrams
      igram(i).t1=igram(i).t1 - firsttime;
      igram(i).t2=igram(i).t2 - firsttime;
       t(i,:)=[igram(i).t1 igram(i).t2];
   end
   datelist=unique(t) ;
   N_timevector=length(datelist);
else
   datelist=[];
   N_timevector=[];
end
   
