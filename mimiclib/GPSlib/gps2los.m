function [GPS]=gps2los(GPS,lookangle)

%
% Compute L.O.S. component of GPS on a structure (see Read_GPS.m for format)
%
% [losrate]=gps2los(GPS)
%
%
% N. Gourmelen, March 2007
%

if nargin == 1
	vect = [cos(11/360*2*pi)*cos(67/360*2*pi) -cos(79/360*2*pi)*cos(67/360*2*pi) cos(23/360*2*pi)];
elseif isstruct(lookangle)
    vect = GenerateLosvector(lookangle(1).sat,lookangle(1).y_first,lookangle(1).sat_height);
else
	vect = lookangle;
end

for ni=1:length(GPS)

    if isfield(GPS(ni),'e_rate')

        GPS(ni).los=[GPS(ni).e_rate GPS(ni).n_rate GPS(ni).height_rate]*vect';GPS(ni).horlos=[GPS(ni).e_rate GPS(ni).n_rate]*vect(1:2)';
        GPS(ni).error_los=sqrt(GPS(ni).e_error.^2*cos(11/360*2*pi).^2*cos(67/360*2*pi).^2+GPS(ni).n_error.^2*cos(79/360*2*pi).^2*cos(67/360*2*pi).^2+GPS(ni).height_error.^2*cos(23/360*2*pi).^2);

    elseif isfield(GPS(ni),'e')

        GPS(ni).los=[GPS(ni).e GPS(ni).n GPS(ni).u]*vect';GPS(ni).horlos=[GPS(ni).e GPS(ni).n]*vect(1:2)';
        if isfield(GPS(ni),'error_e')
            GPS(ni).error_los=sqrt(GPS(ni).error_e.^2*cos(11/360*2*pi).^2*cos(67/360*2*pi).^2+GPS(ni).error_n.^2*cos(79/360*2*pi).^2*cos(67/360*2*pi).^2+GPS(ni).error_u.^2*cos(23/360*2*pi).^2);
        end
        
    end
    
    GPS(ni).lookangle=vect;
    
end


