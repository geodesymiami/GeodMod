function J = date2j(yy,mm,dd)
%DATE2J   J = date2j(yy,mm,dd)
%  
%Returns the Modified Julian Date (MJD) for a given year, month and day.
%If only one input is passed, that input is assumed to be an nx3 matrix,
%where the first, second, and third columns hold the year, month, and day,
%respectively (if fourth, fifth, or sixth columns are present, they are taken
%to be hours, minutes, and seconds).  1900 is added to the year values if
%they are two digits and greater than 80, otherwise 2000 is added.

%
%Note: add 2400000.5 to the MJD to get the Julian Date (JD).

% Do some argument checking
	if nargin == 0
		help date2j
   	return
	end
   
   extra=0;
	if nargin == 1
      if size(yy,2)>3
        extra=yy(:,4)*3600;
        if size(yy,2)>4
         extra=extra+yy(:,5)*60;
           if size(yy,2)>5
         	extra=extra+yy(:,6);
           end
        end
      end
      mm=yy(:,2);
   	dd=yy(:,3);
   	yy=yy(:,1);
	end


	yy=yy+(yy<100)*1900+(yy<80)*100;
   J=floor(365.25*(yy-(mm<=2)))+floor(30.6001*(mm+1+12*(mm<=2)))+dd-679019+extra/86400;
