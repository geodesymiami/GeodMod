function [Stacks]=StackRandFalk(igram,varargin);
%StackRandFalk          - random stacks (modified from Noel's Stacks)
%
%  usage:  [Stacks]=Stack(igram,'NbStck',10,'rating',9999,'repeatNb',1,'LowTmpThresh',0,'HighTmpThresh',9999,'CSpaceLim',[100 200 100 200])
%          plot_igram(igram(2:3))
%
%          Input:   igram:   1xN structure array containing N interferograms with:
%
%          optional arguments (given as 'CLim',[-1,1],'LocInd',[2,2],...)
%
%          'NbStck'        Number of random stacks                   [default   10]
%          'rating'        Maximum rate allow                        [default 9999]
%          'repeatNb'      How many times the same image             [default    1]
%          'MinNbInt'      Minimum number of Interferograms in stack [default    2]
%          'LowTmpThresh'  Minimum time span between images(in yr)   [default    0]
%          'HighTmpThresh' Maximum time span between images(in yr)   [default 9999]
%          'CSpaceLim'     Window for calculating colorscale [x y]   [default  all]          
%          'CLim'          Colorscale bounds                         [default min/max]
%          'StackSortMethod'Method to sort the stack                 [default 'MeanRating']
%                           (MeanRating,TotalTime)                   
%          'DoPlotIgram'     Plot Interferograms ('on','off')          [default on]
%          'PlotStack'     Plot Stacks ('on','off')                  [default on]
%          'MultiFact'      Factor to multiply the stack by           [default 1]
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
    if     strmatch('NbStck',varargin{i})           NbStck=varargin{i+1} ;
    elseif strmatch('rating',varargin{i})           rating=varargin{i+1}  ;
    elseif strmatch('repeatNb',varargin{i})         repeatNb=varargin{i+1}  ;    Coord='off' ;
    elseif strmatch('LowTmpThresh',varargin{i})     LowTmpThresh=varargin{i+1}  ;
    elseif strmatch('HighTmpThresh',varargin{i})    HighTmpThresh=varargin{i+1} ;
    elseif strmatch('CSpaceLim',varargin{i})        CSpaceLim=varargin{i+1} ;
    elseif strmatch('CLim',varargin{i})             CLim=varargin{i+1} ;
    elseif strmatch('DoPlotIgram',varargin{i})      DoPlotIgram=varargin{i+1} ;
    elseif strmatch('PlotStack',varargin{i})        PlotStack=varargin{i+1} ;
    elseif strmatch('StackSortMethod',varargin{i})  StackSortMethod=varargin{i+1} ;
    elseif strmatch('MinNbInt',varargin{i})         MinNbInt=varargin{i+1} ;
    elseif strmatch('Fac',varargin{i})              Fac=varargin{i+1} ;
    elseif strmatch('DirOut',varargin{i})           DirOut=varargin{i+1} ;
    elseif strmatch('Rate',varargin{i})             Rate=varargin{i+1} ;
    else
        error('unknown argument: %s',varargin{i});
    end
end
end

if exist('NbStck')==0 NbStck=10; end
if exist('rating')==0 rating=9999; end
if exist('repeatNb')==0 repeatNb=1; end
if exist('LowTmpThresh')==0 LowTmpThresh=0; end
if exist('HighTmpThresh')==0 HighTmpThresh=9999; end
if exist('DoPlotIgram')==0 DoPlotIgram='off'; end
if exist('PlotStack')==0 PlotStack='off'; end
if exist('MinNbInt')==0 MinNbInt=1; end
if exist('Fac')==0 Fac=1; end % see above MultiFact?
if exist('StackSortMethod')==0 StackSortMethod='MeanRating'; end
if exist('CSpaceLim')==0
  CSpaceLim=[1 size(igram(1).data,2) 1 size(igram(1).data,1)]; 
end

GoodIgramListSaved=zeros(5000,NbStck);
L=length(igram);
k=1;
for iii=1:NbStck
  logmessage(sprintf('Stack Number: %d', iii))
  StackInd=randint(L,1,L,iii)+1;
  StackInd=indivval(StackInd);
  LL=size(StackInd,1);
  total_time_span=0;
  N=0;
  INDICES=[];
  INT=zeros(size(igram(1).data));
  IgramList=[];


%Get date list
%%%%%%%%%%%%%%
DATES=[];rmind=[];
  for I= 1:LL

    IndI=StackInd(I);

    if isfield(igram,'rating')
      Rates=igram(IndI).rating;   
    else
      Rates=1;
    end

    if Rates <= rating

      DATES=[DATES;igram(IndI).t1;igram(IndI).t2];
      
    else

      rmind=[rmind;I];

    end

  end

  intddd=[];

  if exist('Im2Rm')

    for I=1:size(Im2Rm,1)

      ddd=find(Im2Rm(I)==DATES);

      for pl=1:size(ddd,1)
 
        intddd=[intddd;ceil(ddd(pl)/2)];

        if isodd(ddd(pl))
        
          DATES(ddd(pl)+1)=[];DATES(ddd(pl))=[];

        else

          DATES(ddd(pl))=[];DATES(ddd(pl)-1)=[];

        end

      end

    end

  end

  rmind=[rmind;intddd];

  StackInd(rmind)=[];

%Get list of individual dates and shuffle it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  indDATES=indivval(DATES);
  rdind=randperm(size(indDATES,1));


%For each date, check wether it appears more than threshold and sort depending of baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UsedDATES=zeros(size(indDATES));
badStackInd=[];

for DD=1:size(StackInd,1)

    Date1=igram(StackInd(DD)).t1;
    Date2=igram(StackInd(DD)).t2;

    Date1Ind=find(indDATES==Date1);
    Date2Ind=find(indDATES==Date2);

    UsedDATES(Date1Ind)=UsedDATES(Date1Ind)+1;
    UsedDATES(Date2Ind)=UsedDATES(Date2Ind)+1;

    if UsedDATES(Date1Ind) > repeatNb | UsedDATES(Date2Ind) > repeatNb

        badStackInd=[badStackInd;DD];
		repeatNb;

	end
end

IgramList=StackInd;
IgramList(badStackInd)=[];


%DATES
  GoodIgramList=[];

  for StckInd=1:size(IgramList,1)

    t1=igram(IgramList(StckInd)).t1;
    t2=igram(IgramList(StckInd)).t2;

    igram_time_span=(t2-t1)/365.25;

    if isfield(igram,'tforstack')
      igram_time_span=igram(IgramList(StckInd)).tforstack;
    end

    if LowTmpThresh < abs(igram_time_span) & abs(igram_time_span) < HighTmpThresh
 
      GoodIgramList=[GoodIgramList;IgramList(StckInd)];
      if isfield(igram,'rating')
        InterfList(N+1)=str2list(['igram-nb:',num2str(IgramList(StckInd)),'---Rating:',num2str(igram(IgramList(StckInd)).rating),'--',num2str(igram(IgramList(StckInd)).date1),'-',num2str(igram(IgramList(StckInd)).date2)]);
      else
        InterfList(N+1)=str2list(['igram-nb:',num2str(IgramList(StckInd)),'---Rating:',num2str(2),'--',num2str(igram(IgramList(StckInd)).date1),'-',num2str(igram(IgramList(StckInd)).date2)]);
      end
      Int=igram(IgramList(StckInd)).data;

      if igram_time_span > 0;

        INT=INT+Int;
        N=N+1;
      else
  
        INT=INT-Int;
        N=N+1;
      end
    
      total_time_span=total_time_span+abs(igram_time_span);
      INDICES=[INDICES;IndI];
      
      clear Int
      clear date1 ind1 date2 

    else 
    end
  end 

  redond=0;

  for GG=1:iii
	TmpGoodIgramList=zeros(5000,1);TmpGoodIgramList(1:size(GoodIgramList,1),1)=GoodIgramList;
	Tmpp=sum(GoodIgramListSaved(:,GG)-TmpGoodIgramList);
	if Tmpp==0
		redond=1;
    end
  end
 if redond==1
	redond ;
	GoodIgramListSaved(1:size(GoodIgramList,1),iii)=sort(GoodIgramList);

 else

  GoodIgramListSaved(1:size(GoodIgramList,1),iii)=sort(GoodIgramList);


  if size(INDICES,1) < MinNbInt
  else


%INT(240:320,180:220)=NaN;INT(245:310,270:320)=NaN;INT(355:380,235:265)=NaN;

% if N>2
  
%  test=sort(GoodIgramList);
  
%  if k>1
 
%    for ik=1:k

%      if (size(test,1)==size(S(ik).IgramIndices,1))
%      testM=[S(ik).IgramIndices    
        
%      else
      
%      end

%    end

    
    %INT=rmplane(INT,1);
 
    Years=[];
    for IndI=1:size(GoodIgramList,1)
      Years=[Years;str2num(igram(GoodIgramList(IndI)).date1);str2num(igram(GoodIgramList(IndI)).date2)]; 
    end
    indYears=indivval(Years);
    
    if Rate
       INT=INT/total_time_span*Fac ;       % divide by time span to convert from displacment into velocity (rate)
    else
       INT=INT*Fac/length(GoodIgramList);  % keep displacement
    end
    
    Stacks(k).data         = INT;
    Stacks(k).IgramIndices = GoodIgramList;
    Stacks(k).date1        = num2str(indYears(1));
    Stacks(k).date2        = num2str(indYears(size(indYears,1)));
    Stacks(k).x_first      = igram(1).x_first;
    Stacks(k).y_first      = igram(1).y_first;
    Stacks(k).x_step       = igram(1).x_step;
    Stacks(k).y_step       = igram(1).y_step;
    Stacks(k).TotalTime    = total_time_span;
    Stacks(k).x_unit       = igram(1).x_unit;

    k=k+1;

  clear IgramList;
  end
end
end

% Sort Stacks
%%%%%%%%%%%%%
if length(Stacks)==0 
   error ('no data in stack. Try to change MinNbInt to 1 -- existing')
end

% Sort Stacks
%%%%%%%%%%%%%

for i=1:length(Stacks)
	Stacks(i).IgramIndices=sort(Stacks(i).IgramIndices);
	Sizes(i)=size(Stacks(i).IgramIndices,1);
	for j=1:Sizes(i)
        if isfield(igram,'rating')
	  	   tmp(j)=igram(Stacks(i).IgramIndices(j)).rating;
        else
           tmp(j)=1 ;                                              % if no rating given all igrams are assigned rating=1 so that sorting works
        end                                                        % this is not very good.I'd rather should do certain sorting only if rating is given
	end
	Stacks(i).MeanRating=mean(tmp);
end

SizesI=fliplr(indivval(Sizes));

SortedStacks=Stacks;
strind=1;

for i=1:size(SizesI,2)
	clear rind
	ind=find(Sizes==SizesI(i));
	for h=1:size(ind,2)
		rind(h)=Stacks(ind(h)).MeanRating;
	end
	[sortedr,isortedr]=sort(rind);
	SortedStacks(strind:strind+size(ind,2)-1)=Stacks(ind(isortedr));
	strind=strind+size(ind,2);
end
	
Stacks=SortedStacks;clear SortedStacks
%
% Falk's change: Add string with IgramDates 
%

for i=1:length(Stacks)
    for j=1:length(Stacks(i).IgramIndices)
        Stacks(i).IgramDates{j}=[ igram(Stacks(i).IgramIndices(j)).date1 '-'  igram(Stacks(i).IgramIndices(j)).date2 ];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SORT: Falk's Changes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: Removing unnecessary stacks and sorting should be done prior to stack generation
%       to save cpu time and memory 

if strcmp(StackSortMethod,'MeanRating')

%
%  sort according to MeanRating (in some cases Noels version did not work well)
%
  logmessage(sprintf('StackSortMethod: %s',StackSortMethod));

  clear MeanRating sortedStacks 
  [MeanRating{1:length(Stacks)}] = deal(Stacks.MeanRating) ; MeanRating=cell2mat(MeanRating);
  [junk,sortind]=sort(MeanRating,'ascend');
   sortedStacks=Stacks(sortind);
   Stacks=sortedStacks;
%
%  sort according to TotalTime
%
  clear MeanRating sortedStacks
  [MeanRating{1:length(Stacks)}] = deal(Stacks.MeanRating) ; MeanRating=cell2mat(MeanRating);
  meanratinglist=unique(MeanRating);
  for i=1:length(meanratinglist)
      clear ind tmpStacks sortedtmpStacks TotalTime
      ind=find(MeanRating==meanratinglist(i));
      tmpStacks=Stacks(ind);
      [TotalTime{1:length(tmpStacks)}] = deal(tmpStacks.TotalTime) ; TotalTime=cell2mat(TotalTime);
      [junk,sortind]=sort(TotalTime,'descend');
      sortedtmpStacks=tmpStacks(sortind);
      Stacks(ind)=sortedtmpStacks;
  end

elseif strcmp(StackSortMethod,'TotalTime')

%
% sort according to TotalTime  (TODO: after sorting according to TotalTime should
%                                     sort according to MeanRating similar as above )
%
  logmessage(sprintf('StackSortMethod: %s',StackSortMethod));

  clear TotalTime sortedStacks 
  [TotalTime{1:length(Stacks)}] = deal(Stacks.TotalTime) ; TotalTime=cell2mat(TotalTime);
  [junk,sortind]=sort(TotalTime,'descend');
   sortedStacks=Stacks(sortind);
   Stacks=sortedStacks;

else

  error('StackSortMethod %s not recognized',StackSortMethod)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val]=indivval(vect);

clear val
vect=sort(vect);
j=1;
val=zeros(size(vect));
val(1)=vect(1);

for i=2:max(size(vect))
  if vect(i)~=vect(i-1)
   j=j+1;
   val(j)=vect(i);
  else
  end
end

val=val(1:j);

