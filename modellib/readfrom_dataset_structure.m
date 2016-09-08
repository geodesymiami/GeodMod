%script   readfrom_dataset_structure 
% script to read the dataset structure for plotting with plot_Model

datasetq=dataset;   % there seems to be a bug in matlab versions new than 2006b. 

if isfield(datasetq,'x_posting') && ~isempty(datasetq(1).x_posting)
  x_posting       = datasetq(1).x_posting       ;
  y_posting       = datasetq(1).y_posting       ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% bring data into format for DistOp %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [D_1,D_2,D_3,D_4,D_5,D_6,D_7,D_8]=deal(0)         ;

  i=1 ;
  D_1           = true                     ;
  Ndata_1       = datasetq(i).Ndata       ;
  datavec_1     = datasetq(i).datavec     ;
  coord_1       = datasetq(i).coord       ;
  DataSet_1     = datasetq(i).DataSet     ;
  radarlook_1   = datasetq(i).radarlook   ;
  SAR_1         = datasetq(i).SAR       ;
  GPS_1         = datasetq(i).GPS       ;
  sigphi_1      = datasetq(i).sigphi      ;
  if SAR_1
  cx_1          = datasetq(i).cx          ;
  cy_1          = datasetq(i).cy          ;
  data_1        = datasetq(i).data    ;
  amp_1         = ones(size(datasetq(i).data)); 
  amp_1(isnan(datasetq(i).data))=nan;
  if isfield(datasetq(i),'data_mod') data_mod_1=datasetq(i).data_mod; else data_mod_1=datasetq(i).data    ; end
  end


if length(datasetq)>=2
  i=i+1 ;
  D_2           = true                     ;
  Ndata_2       = datasetq(i).Ndata       ;
  datavec_2     = datasetq(i).datavec ;
  coord_2       = datasetq(i).coord       ;
  DataSet_2     = datasetq(i).DataSet     ;
  radarlook_2   = datasetq(i).radarlook   ;
  sigphi_2      = datasetq(i).sigphi      ;
  SAR_2         = datasetq(i).SAR       ;
  GPS_2         = datasetq(i).GPS       ;
  if SAR_2
  cx_2          = datasetq(i).cx          ;
  cy_2          = datasetq(i).cy          ;
  data_2        = datasetq(i).data    ;
  amp_2         = ones(size(datasetq(i).data)); 
  amp_2(isnan(datasetq(i).data))=nan;
  if isfield(datasetq(i),'data_mod') data_mod_2=datasetq(i).data_mod; else data_mod_2=datasetq(i).data    ; end
  end
end

if length(datasetq)>=3
  i=i+1 ;
  D_3           = true                     ;
  Ndata_3       = datasetq(i).Ndata       ;
  datavec_3     = datasetq(i).datavec ;
  coord_3       = datasetq(i).coord       ;
  DataSet_3     = datasetq(i).DataSet     ;
  radarlook_3   = datasetq(i).radarlook   ;
  sigphi_3      = datasetq(i).sigphi      ;
  SAR_3         = datasetq(i).SAR       ;
  GPS_3         = datasetq(i).GPS       ;
  if SAR_3
  cx_3          = datasetq(i).cx          ;
  cy_3          = datasetq(i).cy          ;
  data_3        = datasetq(i).data    ;
  amp_3         = ones(size(datasetq(i).data)); 
  amp_3(isnan(datasetq(i).data))=nan;
  if isfield(datasetq(i),'data_mod') data_mod_3=datasetq(i).data_mod; else data_mod_3=datasetq(i).data    ; end
  end
end

if length(datasetq)>=4
  i=i+1 ;
  D_4           = true                    ;
  Ndata_4       = datasetq(i).Ndata       ;
  datavec_4     = datasetq(i).datavec     ;
  coord_4       = datasetq(i).coord       ;
  DataSet_4     = datasetq(i).DataSet     ;
  radarlook_4   = datasetq(i).radarlook   ;
  cx_4          = datasetq(i).cx          ;
  cy_4          = datasetq(i).cy          ;
  data_4        = datasetq(i).data    ;
  sigphi_4      = datasetq(i).sigphi      ;
  SAR_4         = datasetq(i).SAR       ;
  GPS_4         = datasetq(i).GPS       ;
  amp_4         = ones(size(datasetq(i).data)); 
  amp_4(isnan(datasetq(i).data))=nan;
  if isfield(datasetq(i),'data_mod') data_mod_4=datasetq(i).data_mod; else data_mod_4=datasetq(i).data    ; end
end

if length(datasetq)>=5
  i=i+1 ;
  D_5           = true                    ;
  Ndata_5       = datasetq(i).Ndata       ;
  datavec_5     = datasetq(i).datavec     ;
  coord_5       = datasetq(i).coord       ;
  DataSet_5     = datasetq(i).DataSet     ;
  radarlook_5   = datasetq(i).radarlook   ;
  cx_5          = datasetq(i).cx          ;
  cy_5          = datasetq(i).cy          ;
  sigphi_5      = datasetq(i).sigphi      ;
  SAR_5         = datasetq(i).SAR       ;
  GPS_5         = datasetq(i).GPS       ;
  amp_5         = ones(size(datasetq(i).data)); 
  amp_5(isnan(datasetq(i).data))=nan;
  if isfield(datasetq(i),'data_mod') data_mod_5=datasetq(i).data_mod; else data_mod_5=datasetq(i).data    ; end
end

if length(datasetq)>=6
  i=i+1 ;
  D_6           = true                     ;
  Ndata_6       = datasetq(i).Ndata       ;
  datavec_6     = datasetq(i).datavec     ;
  coord_6       = datasetq(i).coord       ;
  DataSet_6     = datasetq(i).DataSet     ;
  radarlook_6   = datasetq(i).radarlook   ;
  cx_6          = datasetq(i).cx          ;
  cy_6          = datasetq(i).cy          ;
  data_6        = datasetq(i).data    ;
  sigphi_6      = datasetq(i).sigphi      ;
  SAR_6         = datasetq(i).SAR       ;
  GPS_6         = datasetq(i).GPS       ;
  amp_6         = ones(size(datasetq(i).data)); 
  amp_6(isnan(datasetq(i).data))=nan;
  if isfield(datasetq(i),'data_mod') data_mod_6=datasetq(i).data_mod; else data_mod_6=datasetq(i).data    ; end
end

if length(datasetq)>=7
  i=i+1 ;
  D_7           = true                     ;
  Ndata_7       = datasetq(i).Ndata       ;
  datavec_7     = datasetq(i).datavec     ;
  coord_7       = datasetq(i).coord       ;
  DataSet_7     = datasetq(i).DataSet     ;
  radarlook_7   = datasetq(i).radarlook   ;
  cx_7          = datasetq(i).cx          ;
  cy_7          = datasetq(i).cy          ;
  data_7        = datasetq(i).data    ;
  sigphi_7      = datasetq(i).sigphi      ;
  SAR_7         = datasetq(i).SAR       ;
  GPS_7         = datasetq(i).GPS       ;
  amp_7         = ones(size(datasetq(i).data)); 
  amp_7(isnan(datasetq(i).data))=nan;
  if isfield(datasetq(i),'data_mod') data_mod_7=datasetq(i).data_mod; else data_mod_7=datasetq(i).data    ; end
end

if length(datasetq)>=8
  i=i+1 ;
  D_8           = true                     ;
  Ndata_8       = datasetq(i).Ndata       ;
  datavec_8     = datasetq(i).datavec     ;
  coord_8       = datasetq(i).coord       ;
  DataSet_8     = datasetq(i).DataSet     ;
  radarlook_8   = datasetq(i).radarlook   ;
  cx_8          = datasetq(i).cx          ;
  cy_8          = datasetq(i).cy          ;
  data_8        = datasetq(i).data    ;
  sigphi_8      = datasetq(i).sigphi      ;
  SAR_8         = datasetq(i).SAR       ;
  GPS_8         = datasetq(i).GPS       ;
  amp_8         = ones(size(datasetq(i).data)); 
  amp_8(isnan(datasetq(i).data))=nan;
  if isfield(datasetq(i),'data_mod') data_mod_8=datasetq(i).data_mod; else data_mod_8=datasetq(i).data    ; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% bring data into format for anneal %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FA 7/2008: below does the same thing but is much easier
%tmp=cell(length(datasetq),1);;
%[tmp{:}]=deal(datasetq.datavec);  d=cell2mat(tmp');
%[tmp{:}]=deal(datasetq.coord);        coord=cell2mat(tmp');
%[tmp{:}]=deal(datasetq.radarlook);    radarlook=cell2mat(tmp');  
%[tmp{:}]=deal(datasetq.Ndata);        datind=cell2mat(tmp);      datind=cumsum(datind);
%[tmp{:}]=deal(datasetq.DataSet);      DataSet=tmp; 

   d         = [datasetq.datavec]';
   coord     = [datasetq.coord];
   radarlook = [datasetq.radarlook];
   Ndata     = [datasetq.Ndata]';
   sigphi    = [datasetq.sigphi]';
   DataSet   = {datasetq.DataSet}';
   datind    = cumsum(Ndata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% see what GPS are given and put data back into GPS structure %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GPShorz=~isempty(strmatch('GPSeast',DataSet)) & ~isempty(strmatch('GPSnorth',DataSet)) ;
GPSvert=~isempty(strmatch('GPSvert',DataSet)) ;

%ieast =strmatch('GPSeast' ,DataSet);
%inorth=strmatch('GPSnorth',DataSet);
%iup   =strmatch('GPSup',DataSet);
%GPS.xy=datasetq(ieast).coord;
%GPS.enu=zeros(3,length(GPS.xy)); 
%GPS.sig=zeros(3,length(GPS.xy)); 
%if ieast  GPS.enu(1,:)=datasetq(ieast ).datavec; GPS.sig(1,:)=datasetq(ieast ).normalization; end
%if inorth GPS.enu(2,:)=datasetq(inorth).datavec; GPS.sig(2,:)=datasetq(inorth).normalization; end
%if iup    GPS.enu(3,:)=datasetq(iup   ).datavec; GPS.sig(3,:)=datasetq(iup   ).normalization; end
