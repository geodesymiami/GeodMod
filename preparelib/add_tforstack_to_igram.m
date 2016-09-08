function [igram]=add_tforstack_to_igram(igram,date)
%add_tforstack_to_igram - adds tforstack to igram structure
% usage:  [igram]=add_tforstack_to_igram(igram,'date')
%         [igram]=add_tforstack_to_igram(igram,'20020501')
%          Input:   igram:          1xN structure array containing N interferograms with:
%
%         Output:   igram.tforstack   time in yrs since 'date'
%
%  Part of the TimeSeries suite
%  FA, March 2005,   
%  FA, Sep 2006, fixed for cases where mjdstart<mjd1   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ATTENTION date='20020501':'20041201'  not yet incorporated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datestart=date;
mjdstart=date2j(str2num(datestart(3:4)),str2num(datestart(5:6)),str2num(datestart(7:8)));

for j=1:length(igram)
   date1=igram(j).date1;
   date2=igram(j).date2;
   mjd1=date2j(str2num(date1(3:4)),str2num(date1(5:6)),str2num(date1(7:8)));
   mjd2=date2j(str2num(date2(3:4)),str2num(date2(5:6)),str2num(date2(7:8)));
   if mjdstart<=mjd1
      tforstack=(mjd2-mjd1)/365.25  ;
   else 
      tforstack=(mjd2-mjdstart)/365.25  ;
   end
   if tforstack < 0 tforstack=0; end
   igram(j).tforstack=tforstack;
end
return
