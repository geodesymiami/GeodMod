function [igram,n_igrams]=RemoveBadIgrams(igram,maxthresh)
%RemoveBadIgrams        - remove poorly unwrapped interferograms from igram structure
%
% remove interferograms with too much non-unwrapped points
%       maxthresh : remove if non-unwrapped percentage >= thresh
%
dim_data=length(igram) ;

ib=1; ind_badones=zeros(dim_data,1);

for i=1:dim_data
    i_nan=find(isnan(igram(i).data));
    perc=length(i_nan)/prod(size(igram(1).data)) ;
    if (perc >= maxthresh) ind_badones(ib)=i; ib=ib+1; end

   str=sprintf('File %d:  percentage not-unwrapped: %4.2f',i,perc) ; 
   if isfield(igram(i),'date1')
      str=[ str ': ' igram(i).date1(3:end) '-' igram(i).date2(3:end) ',' ];
   end
   logmessage(str)
end

ind_badones(find(ind_badones==0))=[];
igram(ind_badones)=[];
n_igrams=length(igram);
if n_igrams==0
  error('user error: no interferograms left after applying coherence threshold -- exiting') ;
end
