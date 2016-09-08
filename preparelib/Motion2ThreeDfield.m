function []=Motion2ThreeDField(opt)
%Rates2ThreeDField      - Inverts multiple motion files into 3-D velocity field
%
%  usage:  Motion2ThreeDField(motion_list,out_name)
%
%          Options given as opt.dir_in, etc.
%
%          'dir_in'            Input directory (contains e.g. RsatA3)       [default 'off']
%          'out_name'          name to save file                            
%          'igram_name'        filename for igram file                      [default is obtained with extract_name_fromSOdir]
%          'plotdataopt'      plot options                                 [default 'off' ]
%
%  changed to accept another mask.
%
%  V1.0  Falk Amelung, September 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultopt=struct(                                                         ...
        'DoIt'               ,        'on'                    ,            ...
        'in_list'            ,        'off'                   ,            ...
        'out_name'           ,        'off'                   ,            ...
        'NorthComponent'     ,        'off'                   ,            ...
        'plotdataopt'       ,        'off'      )            ;
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],defaultopt); end ;
[opt]=process_defaultoptions(opt,defaultopt);  %display(opt)
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end
global dir_out;
out_name=fullfile(dir_out,'motion_threedfield');
if  ~DoIt                         return; end
if  CheckInOut(in_list,out_name)  return; end

if length(in_list)<=1
   logmessage('One or zero data set -- do not make 3D field')
   return
end

%
%Load all the motions into one multi-dimension motion structure
%
for i=1:length(in_list)
    load(in_list{i});
    tmpmotion(i)  =motion;
    totaltime(i)=motion.TotalTime;
end
motion=tmpmotion;

% East-North-Up system
G=[];d=[];
for i=1:length(motion)
    motion(i).radarlook=extract_hardwired_satparameters(motion(i).DataSet,'LOSvector');    % FA Feb2008 after removing radarlook generation from Roipac2Igram
    G = [G;motion(i).radarlook];
    d = [d;motion(i).data(:)'];
end

if ~NorthComponent
   G(:,2)=0;    % assume u_north=0
end

m=pinv(G)*d;

[enu(1),enu(2),enu(3)] = deal(motion(1));

if (isfield(enu,'IgramIndices') && isfield(enu,'TotalTime') && isfield(enu,'IgramDates'))
    enu            = rmfield(enu,{'IgramIndices' 'TotalTime' 'IgramDates'});
else
    enu            = rmfield(enu,{'TotalTime'});
end

enu(1).DataSet     = 'East' ;
enu(2).DataSet     = 'North';
enu(3).DataSet     = 'Up'   ;

enu(1).data(:)     = m(1,:) ;
enu(2).data(:)     = m(2,:) ;
enu(3).data(:)     = m(3,:) ;

[enu(:).TotalTime] = deal(sum(totaltime));   % summed total time

tmp_enu            = enu;
if ~NorthComponent
    tmp_enu(2)=[];           
end

plotdataopt.googleearthopt.DoIt = 'on' ;

logplot('PlotData_wrapper',out_name,tmp_enu,plotdataopt)
save(out_name,'enu')
