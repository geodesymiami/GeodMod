function []=PlotTimeCoverage(in_list,opt)
%PlotTimeCoverage      -  plots the interferogram time for each data set
%
%  usage:  PlotTimeCoverage(list,opt)
%
%          list              list contains names of interferogram motion files
%                            structure needs to contain the field 'IgramDates'
%
%          opt               options structure (could be plottimecoverageopt, not yet implemented except 'DoIt')
%                            could also be plotdataopt
%          'plotdataopt'      plot options                                 [default 'off' ]
%
%  changed to accept another mask.
%
%  V1.0  Falk Amelung,  December 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultopt=struct(                                                         ...
        'DoIt'               ,        'on'                    ,            ...
        'plotdataopt'        ,        'off'      )            ;
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],defaultopt); end ;
[opt]=process_defaultoptions(opt,defaultopt);  %display(opt)
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end
if  ~DoIt                         return; end
if  CheckInOut(in_list,'')  return; end

%
%Load all the motions into one multi-dimension motion structure
%
for i=1:length(in_list)
    load(in_list{i});
    tmpmotion(i)=motion;
end
motion  = tmpmotion; clear tmpmotion
if  ~isfield(motion,'IgramDates')  return; end

% calcualte limits for plot (rounded first and last image)
tmpdates = char([motion(:).IgramDates]);
plot_date_start    = min ( str2num(tmpdates(:,1:4))  )     ;
plot_date_end      = max ( str2num(tmpdates(:,10:13))) + 1 ;

% calculate locations for lines
y=0.4; date1=[];date2=[];
for i=1:length(motion)
    tmpdates = char([motion(i).IgramDates])  ;
    
    tmpdate1 = YYYYmmDD2DecimalYear(tmpdates(:, 1:8));
    tmpdate2 = YYYYmmDD2DecimalYear(tmpdates(:,10:17));
     
    date1    = [date1 ; tmpdate1];    
    date2    = [date2 ; tmpdate2];
    
    tmpy       = -[1:length(tmpdate1)]/5 + min(y) - 0.4;
    ytext(i)   = tmpy(1);
    datetext(i)= tmpdate1(1);
    y        = [y ; tmpy(:)];
end
date=[date1(:) date2(:)];
y(1)=[];

% plot horizontal lines
a=axes('Box','on','XGrid','on','YTickLabel','');
hold on
for i=1:length(y)
    plot(date(i,:),[y(i) y(i)],'LineWidth',4, 'Color',[0 0 0]);
end
xlim([plot_date_start,plot_date_end])
ylim([min(y)-0.2 max(y)+0.2])

%plot DataSet
for i=1:length(motion)
   text(datetext(i),ytext(i),motion(i).DataSet,'HorizontalAlignment','Right','VerticalAlignment','Baseline')
end
hold off

title('Interferogram time coverage')
