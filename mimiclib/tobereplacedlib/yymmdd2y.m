function [Y]=yyddmm2y(date)

ind=str2num(date(1));

if length(date) == 8
    
elseif ind==9
    date=str2num(date);date=num2str(date+19000000);
 
else
    date=str2num(date);date=num2str(date+20000000);

end

Y=str2num(date(1:4))+(str2num(date(5:6))-1)/12+(str2num(date(7:8))-1)/365;
