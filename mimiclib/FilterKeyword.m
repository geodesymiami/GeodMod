function [OutStruct]=FilterKeyword(InStruct,origin,purpose,data)

%
% Function to filter Keyword structure output from ReadKeywordfile into
% suitable structure for using with geodmod, SBAS, ...
%
% [OutStruct]=FilterKeyword(InStruct,origin,purpose,data)
%
% InStruct is output of ReadKeywordfile
%
% origin: 'roi_pac', 'roi_pac_baseline','irea', 'doris', 'gamma', ...
%
% purpose can be 'geodmod', 'SBAS', ...
%
% data: if present, add a datafield to OutStruct
%
% N. Gourmelen - Aug 2007
%

if nargin == 1  origin='roi_pac';  purpose='geodmod';  end
if nargin == 2  purpose='geodmod';  end

if ~isfield(InStruct,'width')

    switch lower(origin)

        case{'roi_pac'}

            switch lower(purpose)

                case{'geodmod'}

                    OutStruct.width=InStruct.WIDTH;
                    OutStruct.file_length=InStruct.FILE_LENGTH;
                    if isfield(InStruct,'X_FIRST')    OutStruct.x_first=InStruct.X_FIRST;   else    OutStruct.x_first=1;                    end
                    if isfield(InStruct,'Y_FIRST')    OutStruct.y_first=InStruct.Y_FIRST;   else    OutStruct.y_first=InStruct.FILE_LENGTH; end
                    if isfield(InStruct,'X_STEP')     OutStruct.x_step=InStruct.X_STEP;     elseif isfield(InStruct,'RANGE_PIXEL_SIZE')
                        OutStruct.x_step=InStruct.RANGE_PIXEL_SIZE;     else   OutStruct.x_step=1;  end
                    if isfield(InStruct,'Y_STEP')     OutStruct.y_step=InStruct.Y_STEP;     elseif isfield(InStruct,'AZIMUTH_PIXEL_SIZE')
                        OutStruct.y_step=InStruct.AZIMUTH_PIXEL_SIZE;   else   OutStruct.y_step=-1; end
                    if isfield(InStruct,'X_UNIT')     OutStruct.x_unit=InStruct.X_UNIT;       else     OutStruct.x_unit='meters';               end
                    
                    d=findstr(InStruct.DATE12,'-');date1=InStruct.DATE12(1:d-1) ;date2=InStruct.DATE12(d+1:end);sd1=num2str(date1) ;  sd2=num2str(date2);
                    mjd1=date2j(str2num(sd1(1:2)),str2num(sd1(3:4)),str2num(sd1(5:6)));mjd2=date2j(str2num(sd2(1:2)),str2num(sd2(3:4)),str2num(sd2(5:6)));
                    if (strcmp(date1(1),'9')==1)   date1=strcat('19',date1) ; else date1=strcat('20',date1); end    %change Jan 2010 to make work with 2010 data
                    if (strcmp(date2(1),'9')==1)   date2=strcat('19',date2) ; else date2=strcat('20',date2); end    %change Jan 2010 to make work with 2010 data

                    OutStruct.date1=date1;OutStruct.date2=date2;OutStruct.t1=mjd1;OutStruct.t2=mjd2;
		
         	    if isfield(InStruct,'HEADING_DEG')  OutStruct.heading = InStruct.HEADING_DEG ;  end
                    if isfield(InStruct,'LAT_REF1') 
                        OutStruct.los.lat_ref1 =InStruct.LAT_REF1; OutStruct.los.lat_ref2 =InStruct.LAT_REF2; OutStruct.los.lat_ref3 =InStruct.LAT_REF3;OutStruct.los.lat_ref4 =InStruct.LAT_REF4;
                    end
                    if isfield(InStruct,'LON_REF1')
                        OutStruct.los.lon_ref1 =InStruct.LON_REF1; OutStruct.los.lon_ref2 =InStruct.LON_REF2; OutStruct.los.lon_ref3 =InStruct.LON_REF3; OutStruct.los.lon_ref4 =InStruct.LON_REF4;
                    end
                    if isfield(InStruct,'LOOK_REF1')
                        OutStruct.los.look_ref1 =InStruct.LOOK_REF1;  OutStruct.los.look_ref2 =InStruct.LOOK_REF2; OutStruct.los.look_ref3 =InStruct.LOOK_REF4;  OutStruct.los.look_ref4 =InStruct.LOOK_REF4;
                    else OutStruct.los='fixed';
                    end
                    
                    if isfield(InStruct,'HEIGHT') OutStruct.sat_height=InStruct.HEIGHT; else OutStruct.sat_height=8e5; end
                    
                case{'sbas'}

                    OutStruct.width=InStruct.WIDTH;
                    OutStruct.file_length=InStruct.FILE_LENGTH;
                    if isfield(InStruct,'X_FIRST')    OutStruct.x_first=InStruct.X_FIRST;   else    OutStruct.x_first=1;                    end
                    if isfield(InStruct,'Y_FIRST')    OutStruct.y_first=InStruct.Y_FIRST;   else    OutStruct.y_first=InStruct.FILE_LENGTH; end
                    if isfield(InStruct,'X_STEP')     OutStruct.x_step=InStruct.X_STEP;     elseif isfield(InStruct,'RANGE_PIXEL_SIZE')
                        OutStruct.x_step=InStruct.RANGE_PIXEL_SIZE;     else   OutStruct.x_step=1;  end
                    if isfield(InStruct,'Y_STEP')     OutStruct.y_step=InStruct.Y_STEP;     elseif isfield(InStruct,'AZIMUTH_PIXEL_SIZE')
                        OutStruct.y_step=InStruct.AZIMUTH_PIXEL_SIZE;   else   OutStruct.y_step=-1; end
                    if isfield(InStruct,'X_UNIT')     OutStruct.x_unit=InStruct.X_UNIT;       else     OutStruct.x_unit='meters';               end

                    if isfield(InStruct,'DATE12')
                        d = findstr(InStruct.DATE12,'-') ;
                        date1 = InStruct.DATE12(1:d-1)   ;
                        date2 = InStruct.DATE12(d+1:end) ;
                        sd1 = num2str(date1) ;  
                        sd2 = num2str(date2) ;
                        mjd1=date2j(str2num(sd1(1:2)),str2num(sd1(3:4)),str2num(sd1(5:6)));mjd2=date2j(str2num(sd2(1:2)),str2num(sd2(3:4)),str2num(sd2(5:6)));
                        OutStruct.date1 = date1 ;  OutStruct.date2=date2;OutStruct.t1=mjd1;OutStruct.t2=mjd2;
                        OutStruct.delt=mjd2-mjd1;  
                    end
                    OutStruct.wavelength  = InStruct.WAVELENGTH ;     OutStruct.sat_height   = InStruct.HEIGHT ;
                    OutStruct.near_range  = InStruct.RGE_REF1*1000 ;  OutStruct.far_range    = InStruct.RGE_REF2*1000;OutStruct.near_LookAng=InStruct.LOOK_REF1 ;
                    OutStruct.far_LookAng = InStruct.LOOK_REF2 ;      OutStruct.earth_radius = InStruct.EARTH_RADIUS ;
                    
                otherwise

                    disp('Unknown purpose')

            end

        case{'roi_pac_baseline'}

            switch lower(purpose)

                case{'geodmod'}
OutStruct.los.ul=InStruct;

                case{'sbas'}

                    if isfield(InStruct,'P_BASELINE_TOP_ODR')
                        OutStruct.bperptop=InStruct.P_BASELINE_TOP_ODR;
                        OutStruct.bperpbot=InStruct.P_BASELINE_BOTTOM_ODR;
                    elseif isfield(InStruct,'P_BASELINE_TOP_HDR')
                        OutStruct.bperptop=InStruct.P_BASELINE_TOP_HDR;
                        OutStruct.bperpbot=InStruct.P_BASELINE_BOTTOM_HDR;
                    else
                        OutStruct.bperptop=0;
                        OutStruct.bperpbot=0;
                    end

                otherwise

                    disp('Unknown purpose')

            end

        case{'doris'}

        case{'irea'}
            
            OutStruct.width       =  InStruct.WIDTH       ;
            OutStruct.file_length =  InStruct.FILE_LENGTH ;
            OutStruct.proj        =  InStruct.proj        ;
            OutStruct.x_first     =  InStruct.X_FIRST     ;  
            OutStruct.y_first     =  InStruct.Y_FIRST + (InStruct.FILE_LENGTH-1)*InStruct.Y_STEP ;
            OutStruct.x_step      =  InStruct.X_STEP      ; 
            OutStruct.y_step      = -InStruct.Y_STEP      ;
            OutStruct.x_unit      =  InStruct.X_UNIT      ;
            OutStruct.Unit        =  InStruct.Unit        ;
            if isfield(InStruct,'DATE12')
                        d = findstr(InStruct.DATE12,'-') ;
                        date1 = InStruct.DATE12(1:d-1)   ;
                        date2 = InStruct.DATE12(d+1:end) ;
                        sd1   = num2str(date1) ;  
                        sd2   = num2str(date2) ;
                        mjd1  = date2j(str2num(sd1(1:2)),str2num(sd1(3:4)),str2num(sd1(5:6))) ;
                        mjd2  = date2j(str2num(sd2(1:2)),str2num(sd2(3:4)),str2num(sd2(5:6))) ;
                        OutStruct.date1 = date1 ;  OutStruct.date2 = date2 ;
                        OutStruct.t1   = mjd1 ;  OutStruct.t2 = mjd2 ;
                        OutStruct.delt = mjd2-mjd1 ;  
            end
            if isfield(InStruct,'sat_height') ;  OutStruct.sat_height = InStruct.sat_height ;  else OutStruct.sat_height = 7933882 ; end
            if isfield(InStruct,'wavelength') ;  OutStruct.wavelength = InStruct.wavelength ;  else OutStruct.wavelength = 0.05656 ; end

        case{'gamma'}

        case{'dem'}

            OutStruct.width         = InStruct.WIDTH        ;
            OutStruct.file_length   = InStruct.FILE_LENGTH  ;
            OutStruct.x_first       = InStruct.X_FIRST      ;
            OutStruct.y_first       = InStruct.Y_FIRST      ;
            OutStruct.x_step        = InStruct.X_STEP       ;
            OutStruct.y_step        = InStruct.Y_STEP       ;
            OutStruct.x_unit        = InStruct.X_UNIT       ;
            OutStruct.Unit          = 'm'                   ;

        otherwise

            disp('Unknown origin')

    end

else
    if (strcmp(origin,'geodmod') & strcmp(purpose,'roi_pac'))
        OutStruct.WIDTH         = InStruct.width        ;
        OutStruct.FILE_LENGTH   = InStruct.file_length  ;
        OutStruct.X_FIRST       = InStruct.x_first      ;
        OutStruct.Y_FIRST       = InStruct.y_first      ;
        OutStruct.X_STEP        = InStruct.x_step       ;
        OutStruct.Y_STEP        = InStruct.y_step       ;
    else
        OutStruct=InStruct;
    end
end



