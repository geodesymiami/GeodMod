function [Y]=yyyy_ddd2y(date)

for ni=1:length(date)
    ind=num2str(sprintf('%0.3f',date(ni))) ;
    Y(ni) = str2num(ind(1:4))+(str2num(ind(6:8)))/365;
end