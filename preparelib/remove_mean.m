function [igram]=remove_mean(igram)
%remove_mean            - remove igrams based on given acquisition date
%
%  usage:  [igram]=remove_mean(igram)
%          [igram]=remove_mean(igram)
%
%          removes interferograms with matching date from igram structure
%
%  Part of the TimeSeries suite
%  FA, March 2005,   

for j=1:length(igram)
    tmpdata=igram(j).data;
    tmean=nanmean(nanmean(tmpdata));
    igram(j).data=igram(j).data-tmean;
end

return
