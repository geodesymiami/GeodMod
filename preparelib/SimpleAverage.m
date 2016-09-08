function [Stacks]=SimpleAverage(igram,varargin);
%SimpleAverage          - simple average of multiple interferograms (to be incorporated into Noel's MakeStack.m)
%
%
%  usage:  [Stacks]=StackIntRand3(igram,'Fac',1,'CSpaceLim',[100 200 100 200])
%          plot_igram(igram(2:3))
%
%          Input:   igram:   1xN structure array containing N interferograms with:
%
%          optional arguments (given as 'CLim',[-1,1],'LocInd',[2,2],...)
%
%          'Factor'                                [default 1 ]
%          'rating'        Maximum rate allow                        [default 9999]
%          'repeatNb'      How many times the same image             [default    1]
%          'MinNbInt'      Minimum number of Interferograms in stack [default    2]
%          'LowTmpThresh'  Minimum time span between images(in yr)   [default    0]
%          'HighTmpThresh' Maximum time span between images(in yr)   [default 9999]
%          'CSpaceLim'     Window for calculating colorscale [x y]   [default  all]          
%          'CLim'          Colorscale bounds                         [default min/max]
%          'PlotIgram'     Plot Interferograms ('on','off')          [default on]
%          'PlotStack'     Plot Stacks ('on','off')                  [default on]
%          'Rate'          calculate rates (displacement) (on,off)   [default on]
% 
%  Part of the TimeSeries suite
%  NG, March 2005,
%
% process the arguments
%
disp(['working... ' mfilename ])
if length(varargin) >= 1
if (floor(length(varargin)/2)~=length(varargin)/2) , error ('argument missing') , end
for i=1:2:length(varargin)
    if     strmatch('NbStck',varargin{i})          NbStck=varargin{i+1} ;
    elseif strmatch('Fac',varargin{i})             Fac=varargin{i+1}  ;
    elseif strmatch('rating',varargin{i})       rating=varargin{i+1}  ;
    elseif strmatch('repeatNb',varargin{i})       repeatNb=varargin{i+1}  ;    Coord='off' ;
    elseif strmatch('LowTmpThresh',varargin{i})       LowTmpThresh=varargin{i+1}  ;
    elseif strmatch('HighTmpThresh',varargin{i})      HighTmpThresh=varargin{i+1} ;
    elseif strmatch('CSpaceLim',varargin{i})      CSpaceLim=varargin{i+1} ;
    elseif strmatch('CLim',varargin{i})      CLim=varargin{i+1} ;
    elseif strmatch('PlotIgram',varargin{i})           PlotIgram=varargin{i+1} ;
    elseif strmatch('PlotStack',varargin{i})           PlotStack=varargin{i+1} ;
    elseif strmatch('MinNbInt',varargin{i})      MinNbInt=varargin{i+1} ;
    elseif strmatch('Rate',varargin{i})          Rate=varargin{i+1} ;
    else
        error('unknown argument: %s',varargin{i});
    end
end
end

if exist('Fac')==0           Fac=1; end
if exist('rating')==0        rating=9999; end
if exist('repeatNb')==0      repeatNb=1; end
if exist('LowTmpThresh')==0  LowTmpThresh=0; end
if exist('HighTmpThresh')==0 HighTmpThresh=9999; end
if exist('PlotIgram')==0     PlotIgram='on'; end
if exist('PlotStack')==0     PlotStack='on'; end
if exist('MinNbInt')==0      MinNbInt=2; end
if exist('Rate')==0          Rate=1; end
if exist('CSpaceLim')==0     CSpaceLim=[1 size(igram(1).data,2) 1 size(igram(1).data,1)]; end

int=zeros(size(igram(1).data));TTspan=0;
for i=1:length(igram);

    t1=igram(i).t1;
    t2=igram(i).t2;

    Tsp=(t2-t1)/365.25;

    if isfield(igram(i),'tforstack')
      Tsp=igram(i).tforstack;
    end

    TTspan=TTspan+abs(Tsp);
    int=int+igram(i).data;
end

if Rate
   int=int/TTspan*Fac ;       % divide by time span to convert from displacment into velocity (rate)
else
   int=int*Fac/length(igram) ;              % keep displacement
end

Stacks.data=int;
Stacks.x_first=igram(1).x_first;
Stacks.y_first=igram(1).y_first;
Stacks.x_step =igram(1).x_step;
Stacks.y_step =igram(1).y_step;
Stacks.x_unit =igram(1).x_unit;
Stacks.IgramIndices=[1:length(igram)];
Stacks.TotalTime=TTspan ;

%                           
% Falk's change: Add string with IgramDates 
%          
for j=1:length(Stacks.IgramIndices)
    Stacks.IgramDates{j}=[ igram(Stacks.IgramIndices(j)).date1 '-'  igram(Stacks.IgramIndices(j)).date2 ];
end


