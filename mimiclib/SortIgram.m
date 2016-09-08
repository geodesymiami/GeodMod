function [qigram,n_igrams]=SortIgram(igram)
% SortIgram            -  Sort interferograms in igram structure
%
% Sort interferograms
%
N=length(igram) ;

% first zero out all fields
t=zeros(N,2);
if isfield(igram(1),'t1')
    for i=1:N
        t(i,:)=[igram(i).t1 igram(i).t2];
    end
elseif isfield(igram(1),'t')
    for i=1:N
        t(i,:)=[igram(i).t];
    end
end

[nt,ind]=sortrows(t);

for i=1:N
     nigram(i)=igram(ind(i));                     % changed March 2005
%    nigram(i).data=igram(ind(i)).data;
%    nigram(i).t1=igram(ind(i)).t1;
%    nigram(i).t2=igram(ind(i)).t2;
%    nigram(i).delt=igram(ind(i)).delt;
%    if isfield(igram,'date1') nigram(i).date1=igram(ind(i)).date1; end
%    if isfield(igram,'date2') nigram(i).date2=igram(ind(i)).date2; end
%    if isfield(igram,'par') nigram(i).par=igram(ind(i)).par; end
%    if isfield(igram,'atmo') nigram(i).atmo=igram(ind(i)).atmo; end
%    if isfield(igram,'atmo') nigram(i).atmo=igram(ind(i)).atmo; end
%    if isfield(igram,'noatmo') nigram(i).noatmo=igram(ind(i)).noatmo; end
end

% remove double data and put data back into igram
[junk,ind]=unique(nt,'rows');
for i=1:length(ind)
   qigram(i)=nigram(ind(i));
end
n_igrams=length(qigram);
