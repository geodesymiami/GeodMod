function  [Fac]=convert_unit(out_unit,in_unit,wavelength)

%
%  Output a unit conversion factor
%
%  usage:  [factor]=convert_unit(out_unit,in_unit,wavelength)
%
%       in_unit can be a structure with data and keywords. If so it checks
%       wether fields wavelength and Unit are present. If not it assumes C
%       band and radian unit.
%
% 
% November 2007 Noel Gourmelen
%

%%%
%%%

if nargin < 2  in_unit    = 'radian';  end
if nargin < 3  Wavelength = 0.0566  ;  end
    
if isstruct(in_unit)
    if isfield(in_unit(1),'Unit')        Unit=in_unit(1).Unit;                end
    if isfield(in_unit(1),'wavelength')  Wavelength=in_unit(1).wavelength;    end
else  Unit=in_unit;  
end

if isstruct(out_unit)
    if isfield(out_unit(1),'wavelength')  Wavelength=out_unit(1).wavelength;  end
    if isfield(out_unit(1),'Unit')        out_unit=out_unit(1).Unit;          end  
end

switch(Unit)

    case('radian')
        
        if     strcmp(out_unit,'m/yr')   Fac=Wavelength/4/pi;           rate_flag=true;    %factor for conversion into m/yr
        elseif strcmp(out_unit,'cm/yr')  Fac=Wavelength/4/pi*100;       rate_flag=true;    %factor for conversion into cm/yr
        elseif strcmp(out_unit,'mm/yr')  Fac=Wavelength/4/pi*1000;      rate_flag=true;    %factor for conversion into mm/yr
        elseif strcmp(out_unit,'m')      Fac=Wavelength/4/pi;           rate_flag=true;   %factor for conversion into m
        elseif strcmp(out_unit,'radian') Fac=-1;                        rate_flag=false;    %factor for conversion into mm/yr
        elseif strcmp(out_unit,'unity')  Fac=1;                         rate_flag=false;   %factor for conversion into m
        end

    case('m/yr')
        
        if     strcmp(out_unit,'radian')   Fac=4*pi/Wavelength;         rate_flag=true;    %factor for conversion into m/yr
        elseif strcmp(out_unit,'m/yr')     Fac=1;                       rate_flag=true;    %factor for conversion into cm/yr
        elseif strcmp(out_unit,'cm/yr')    Fac=100;                     rate_flag=true;    %factor for conversion into cm/yr
        elseif strcmp(out_unit,'mm/yr')    Fac=1000;                    rate_flag=true;    %factor for conversion into mm/yr
        elseif strcmp(out_unit,'m')        Fac=1;                       rate_flag=false;   %factor for conversion into m
        end

    case('cm/yr')
        
        if     strcmp(out_unit,'radian')   Fac=4*pi/Wavelength/100;     rate_flag=true;    %factor for conversion into m/yr
        elseif strcmp(out_unit,'m/yr')     Fac=0.01 ;                   rate_flag=true;    %factor for conversion into cm/yr
        elseif strcmp(out_unit,'cm/yr')    Fac=1    ;                   rate_flag=true;    %factor for conversion into cm/yr
        elseif strcmp(out_unit,'mm/yr')    Fac=10   ;                   rate_flag=true;    %factor for conversion into mm/yr
        elseif strcmp(out_unit,'m')        Fac=0.01 ;                   rate_flag=false;   %factor for conversion into m
        end

    case('mm/yr')
        
        if     strcmp(out_unit,'radian')   Fac=4*pi/Wavelength/1000;    rate_flag=true;    %factor for conversion into m/yr
        elseif strcmp(out_unit,'cm/yr')    Fac=0.1;                     rate_flag=true;    %factor for conversion into cm/yr
        elseif strcmp(out_unit,'m/yr')     Fac=0.001;                   rate_flag=true;    %factor for conversion into mm/yr
        elseif strcmp(out_unit,'m')        Fac=0.001;                   rate_flag=false;   %factor for conversion into m
        end
        
    case('m')
        
        if     strcmp(out_unit,'radian')   Fac=4*pi/Wavelength;         rate_flag=true;    %factor for conversion into m/yr
        elseif strcmp(out_unit,'cm/yr')    Fac=100;                     rate_flag=true;    %factor for conversion into cm/yr
        elseif strcmp(out_unit,'mm/yr')    Fac=1000;                    rate_flag=true;    %factor for conversion into mm/yr
        elseif strcmp(out_unit,'m')        Fac=1;                       rate_flag=false;   %factor for conversion into m
        elseif strcmp(out_unit,'unity')    Fac=1;                       rate_flag=false;   %factor for conversion into m
        end
        
    case('cm')
        
        if     strcmp(out_unit,'radian')   Fac=4*pi/Wavelength/100;     rate_flag=true;    %factor for conversion into m/yr
        elseif strcmp(out_unit,'cm')       Fac=1;                       rate_flag=true;    %factor for conversion into cm/yr
        elseif strcmp(out_unit,'mm')       Fac=10;                      rate_flag=true;    %factor for conversion into mm/yr
        elseif strcmp(out_unit,'m')        Fac=0.01;                    rate_flag=false;   %factor for conversion into m
        elseif strcmp(out_unit,'unity')    Fac=1;                       rate_flag=false;   %factor for conversion into m
        end
        
    case('generic')
                                           Fac = 1 ;                    rate_flag=false;   %factor for conversion into m
end
