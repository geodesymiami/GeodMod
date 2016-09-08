function [igram]=add_rating_to_igram(igram,rfile)
%add_rating_to_igram    - add ratings to the igram structure
% usage:  add_rating_to_igram(igram,rfile)
%          add_rating_to_igram(igram,'/RAID1/amelung_as_apex/HawA6B400More2Yrs/igram_rating.log')
%
%          Input:   igram:          1xN structure array containing N interferograms with:
%                   rfile:          ASCII file containing the rating for each interferogram
%                                   in the form (produced by select_and_stack.pl)
%
%                                   PO_IFGRAM_HawA6B400More2Yrs_20010526-20040111       2
%
%                                   searches for 'date1-date2' in the dir name
%
%  Part of the TimeSeries suite
%  FA, March 2005,   

%cc;

[dirs,ratings] = textread(rfile,'%s%d');     % reads ratingfile

for j=1:length(igram)
   datestr=strcat(igram(j).date1,'-',igram(j).date2);
   dirmatch=zeros(length(dirs),1);
   for i=1:length(dirs)
      dirmatch(i)=length(strfind(char(dirs(i)),datestr));
   end
   if length(find(dirmatch==1))==0 error('no match for igram %d : %s',j,datestr) ; end
   if length(find(dirmatch==1))> 1 error('more that one match for igram %d %s',j ,datestr) ; end
   igram(j).rating=ratings(find(dirmatch==1)) ;
end

return
