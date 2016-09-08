function [outname]=ChangeDates(Innames)

%
% change dates format
%
month=Innames(3:5);
if (month == 'JAN')
	monthout='01';
elseif (month == 'FEB')
	monthout='02';
elseif (month == 'MAR')
	monthout='03';
elseif (month == 'APR')
	monthout='04';
elseif (month =='MAY')
	monthout='05';
elseif (month == 'JUN')
	monthout='06';
elseif (month == 'JUL')
	monthout='07';
elseif (month =='AUG')
	monthout='08';
elseif (month =='SEP')
	monthout='09';
elseif (month =='OCT')
	monthout='10';
elseif (month =='NOV')
	monthout='11';
elseif (month=='DEC')
	monthout='12';
end

outname=[Innames(1:2),monthout,Innames(6:7)];
