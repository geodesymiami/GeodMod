function [GPS]=GPStmsr2GPSrate(GPS,opt)

%
% Compute rates from structure of GPS timeseries obtained from Read_GPS.m
%
% [InSAR]=GPStmsr2GPS(opt)
%
% 'GPS'    : GPS structure    
%
% opt can be 
%
%
%             - model      : Model for rate computation. Can be linear (default), annual (Annual + Semi annual), ...
%
%             - time_range : Time range for rate computation. [year_min year_max] 
%			
%
%
% N. Gourmelen - June 2007
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultopt=struct(                                                           ...
		'model'     	                 ,        'linear'      ,            ...
        'time_range'                     ,        'off'         )             ;
    
    
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],[]); end ;
[opt]=process_defaultoptions(opt,defaultopt);  display(opt)
f=fieldnames(opt) ;
for i=1:length(f)
    eval([char(f{i}) '= opt.(f{i}) ;' ]) ;
end

%%%%%%%%%%%%%%%%%%%%%
% Select Time Range %
%%%%%%%%%%%%%%%%%%%%%

if time_range
    for ni = 1:length(GPS)

        tmp_ind = find(GPS(ni).date_years < time_range(1) |  GPS(ni).date_years > time_range(2));


    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read, Clean and Store data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(model)

    case('linear')

        for ni=1:length(GPS)

            tmp_ind = linspace(1,length(GPS(ni).date_years),length(GPS(ni).date_years));
            
            if time_range  tmp_ind = find(GPS(ni).date_years < time_range(1) |  GPS(ni).date_years > time_range(2));  end
            
            if (length(tmp_ind)==length(GPS(ni).date_years))

               disp(sprintf('No GPS measurement in the given time range for station'));

            else

                if isfield(GPS(ni),'e')
                    tmp=pinv([GPS(ni).date_years ones(size(GPS(ni).date_years,1),1)])*GPS(ni).e;
                    GPS(ni).e_rate=tmp(1);
                end
                if isfield(GPS(ni),'n')
                    tmp=pinv([GPS(ni).date_years ones(size(GPS(ni).date_years,1),1)])*GPS(ni).n;
                    GPS(ni).n_rate=tmp(1);
                end

                if isfield(GPS(ni),'u')
                    tmp=pinv([GPS(ni).date_years ones(size(GPS(ni).date_years,1),1)])*GPS(ni).height;
                    GPS(ni).height_rate=tmp(1);
                end
                if isfield(GPS(ni),'los')
                    tmp=pinv([GPS(ni).date_years ones(size(GPS(ni).date_years,1),1)])*GPS(ni).los;
                    GPS(ni).los_rate=tmp(1);
                end
                if isfield(GPS(ni),'horlos')
                    tmp=pinv([GPS(ni).date_years ones(size(GPS(ni).date_years,1),1)])*GPS(ni).los;
                    GPS(ni).horlos_rate=tmp(1);
                end

            end

        end

    case('seasonal')

        for ni=1:length(GPS) % Loop thropugh GPS to get model parameters

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Filter the signal for outliers %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            gpsdates=GPS(ni).date_years';gpsdata=GPS(ni).los';insardates=Datess;insardata=Data(ni,:);

            for ny=1:3 		% Loop 3 times to clean outliers from GPS data
                datamean=mean(gpsdata(isfinite(gpsdata)));
                datastd=std(gpsdata(isfinite(gpsdata)));
                gpsind=find(abs(gpsdata-datamean)<=3*datastd);
                gpsdata=gpsdata(gpsind);gpsdates=gpsdates(gpsind);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Get GPS model parameters %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            GaussianGPS=[ones(size(gpsdates))' gpsdates' sin(2*pi*gpsdates') cos(2*pi*gpsdates') sin(4*pi*gpsdates') cos(4*pi*gpsdates')];
            GaussianInSAR=[ones(size(insardates))' insardates' sin(2*pi*insardates') cos(2*pi*insardates') sin(4*pi*insardates') cos(4*pi*insardates')];

            ModeledInSAR(ni,:)=GaussianInSAR*pinv(GaussianGPS)*gpsdata';

        end

        InSARdiff=Data-ModeledInSAR;
        Gaussian=[InSAR(1).Indices ones(size(Data,1),1)];

        for ny=1:size(Data,2)
            PlaneCoeff(:,ny)=pinv(Gaussian)*InSARdiff(:,ny);
            InSAR(ny).datafit=Data(:,ny)-Gaussian*PlaneCoeff(:,ny);
            InSAR(ny).PlaneCoeff4Fit=PlaneCoeff(:,ny);

            if isstruct(InSARmat)   % Remove the plane to the InSAR time series matrice if given as input

                [XX,YY]=meshgrid(1:size(InSARmat(ny).data,2),1:size(InSARmat(ny).data,1));
                GaussianM=[YY(:),XX(:),ones(size(XX,1)*size(XX,2),1)];
                plane=GaussianM*PlaneCoeff(:,ny);plane=reshape(plane,size(XX,1),size(XX,2));
                InSARmat(ny).datafit=InSARmat(ny).data-plane;

            end

        end

end





  
        %                if isfield(GPS(ni),'lon')
        %                        tmp=pinv([GPS(ni).date_years ones(size(GPS(ni).date_years,1),1)])*GPS(ni).lon;
        %                        GPS(ni).lon_rate=tmp(1);
        %                end
        %                if isfield(GPS(ni),'lat')
        %                        tmp=pinv([GPS(ni).date_years ones(size(GPS(ni).date_years,1),1)])*GPS(ni).lat;
        %                        GPS(ni).lat_rate=tmp(1);
        %               end
       