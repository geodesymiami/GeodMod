function [inverseopt]=InitializeInverseopt(dataset,inverseopt,modelopt);
%MakeBounds   - generates bounds, default values if no bounds are given 
%
%usage: [inverseopt]=MakeBounds(dataset,inverseopt,modelbounds);
%
%  if no sources given then N_disloc is set to 1
%
%  inverseopt is needed as input because it contains the bounds if given
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  we work in Cartesian (xy) coordinates                                                    %
%  TODO: Shall we have the ability to work in different coordinate systems, given by an option %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FA, October, 2006 
%FA, modified June 2007 
%FA, modified August 2007 

f=fieldnames(inverseopt); for i=1:length(f) eval([char(f{i}) '= inverseopt.(f{i}) ;' ]); end
f=fieldnames(modelopt);   for i=1:length(f) eval([char(f{i}) '= modelopt.(f{i}) ;']);    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  Make bounds   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% First transform given bounds into xy coordinate if necessary
%
% TODO: this is a dirty fix that works only for one dislocation and one Mogi source. We need something like:
% if N_disloc  disloc_bounds=FillCoordinates(disloc_bounds,'disloc'); end
% if N_mogi    mogi_bounds  =FillCoordinates(mogi_bounds,  'mogi'  ); end

  if exist('disloc_bounds','var') 
     if ~isfield(disloc_bounds,'xy') && isfield(disloc_bounds,'lola')   
        tmpmodelopt=InitializeModelopt(struct('N_disloc',N_disloc),plotdataopt.basemap);
        disloc_bounds.xy(:,1)=modelpar_lola2xy(disloc_bounds.lola(:,1),plotdataopt.basemap,tmpmodelopt,1);
        disloc_bounds.xy(:,2)=modelpar_lola2xy(disloc_bounds.lola(:,2),plotdataopt.basemap,tmpmodelopt,1);
     end
  end

    if exist('fault_bounds','var') 
     if ~isfield(fault_bounds,'xy') && isfield(fault_bounds,'lola')   
        tmpmodelopt=InitializeModelopt(struct('N_fault',N_fault),plotdataopt.basemap);
        fault_bounds.xy(:,1)=modelpar_lola2xy(fault_bounds.lola(:,1),plotdataopt.basemap,tmpmodelopt,1);
        fault_bounds.xy(:,2)=modelpar_lola2xy(fault_bounds.lola(:,2),plotdataopt.basemap,tmpmodelopt,1);
     end
    end
  
  if exist('mogi_bounds','var') 
     if ~isfield(mogi_bounds,'xy') && isfield(mogi_bounds,'lola')       
        tmpmodelopt=InitializeModelopt(struct('N_mogi',N_mogi),plotdataopt.basemap);
        mogi_bounds.xy(:,1)=modelpar_lola2xy(mogi_bounds.lola(:,1),plotdataopt.basemap,tmpmodelopt,1); 
        mogi_bounds.xy(:,2)=modelpar_lola2xy(mogi_bounds.lola(:,2),plotdataopt.basemap,tmpmodelopt,1); 
     end
  end
  if exist('penny_bounds','var') 
     if ~isfield(penny_bounds,'xy') && isfield(penny_bounds,'lola')       
        tmpmodelopt=InitializeModelopt(struct('N_penny',N_penny),plotdataopt.basemap);
        penny_bounds.xy(:,1)=modelpar_lola2xy(penny_bounds.lola(:,1),plotdataopt.basemap,tmpmodelopt,1); 
        penny_bounds.xy(:,2)=modelpar_lola2xy(penny_bounds.lola(:,2),plotdataopt.basemap,tmpmodelopt,1); 
     end
  end
  if exist('mctigue_bounds','var')
     if ~isfield(mctigue_bounds,'xy') && isfield(mctigue_bounds,'lola')       
        tmpmodelopt=InitializeModelopt(struct('N_mctigue',N_mctigue),plotdataopt.basemap);
        mctigue_bounds.xy(:,1)=modelpar_lola2xy(mctigue_bounds.lola(:,1),plotdataopt.basemap,tmpmodelopt,1); 
        mctigue_bounds.xy(:,2)=modelpar_lola2xy(mctigue_bounds.lola(:,2),plotdataopt.basemap,tmpmodelopt,1); 
     end
  end
  if exist('yang_bounds','var') 
     if ~isfield(yang_bounds,'xy') && isfield(yang_bounds,'lola')       
        tmpmodelopt=InitializeModelopt(struct('N_yang',N_yang),plotdataopt.basemap);
        yang_bounds.xy(:,1)=modelpar_lola2xy(yang_bounds.lola(:,1),plotdataopt.basemap,tmpmodelopt,1); 
        yang_bounds.xy(:,2)=modelpar_lola2xy(yang_bounds.lola(:,2),plotdataopt.basemap,tmpmodelopt,1); 
     end
  end
  if exist('squaredisloc_bounds','var') 
     if ~isfield(squaredisloc_bounds,'xy') && isfield(squaredisloc_bounds,'lola')       
        tmpmodelopt=InitializeModelopt(struct('N_squaredisloc',N_squaredisloc),plotdataopt.basemap);
        squaredisloc_bounds.xy(:,1)=modelpar_lola2xy(squaredisloc_bounds.lola(:,1),plotdataopt.basemap,tmpmodelopt,1); 
        squaredisloc_bounds.xy(:,2)=modelpar_lola2xy(squaredisloc_bounds.lola(:,2),plotdataopt.basemap,tmpmodelopt,1); 
     end
  end
  if exist('multidisloc_bounds','var') 
     if ~isfield(multidisloc_bounds,'xy') && isfield(multidisloc_bounds,'lola')       
        tmpmodelopt=InitializeModelopt(struct('N_multidisloc',N_multidisloc),plotdataopt.basemap);
        multidisloc_bounds.xy(:,1)=modelpar_lola2xy(multidisloc_bounds.lola(:,1),plotdataopt.basemap,tmpmodelopt,1); 
        multidisloc_bounds.xy(:,2)=modelpar_lola2xy(multidisloc_bounds.lola(:,2),plotdataopt.basemap,tmpmodelopt,1); 
     else
         multidisloc_bounds.xy = multidisloc_bounds.xy(1:10+length(multidislocopt.ind),:);%ensures correct length (for testinh there may be longer bounds arrays in the *min file)
     end
  end
  
  if exist('peas_bounds','var') 
    % Not needed?
  end
  

[modelopt,default_bounds] = InitializeModelopt(modelopt,plotdataopt.basemap);

bounds=[];
for i=1:double(N_disloc) 
    if exist('disloc_bounds','var') && length(disloc_bounds.xy) >= 10  bounds = [bounds; disloc_bounds.xy(1:10,:)]; disloc_bounds.xy(1:10,:)=[];  
       else                                                            bounds = [bounds; default_bounds.disloc.xy];
   end
end
for i=1:double(N_fault) 
    if exist('fault_bounds','var')  && length(fault_bounds.xy) >= 9     bounds = [bounds; fault_bounds.xy(1:9,:)]; fault_bounds.xy(1:9,:)=[];  
    else                                                                bounds = [bounds; default_bounds.fault.xy];
   end
end
for i=1:double(N_mogi) 
    if exist('mogi_bounds','var') && length(mogi_bounds.xy) >= 4       bounds = [bounds; mogi_bounds.xy(1:4,:)]; mogi_bounds.xy(1:4,:)=[];  
       else                                                            bounds = [bounds; default_bounds.mogi.xy];
   end
end
for i=1:double(N_penny) 
    if exist('penny_bounds','var') && length(penny_bounds.xy) >= 5     bounds = [bounds; penny_bounds.xy(1:5,:)]; penny_bounds.xy(1:5,:)=[];  
       else                                                            bounds = [bounds; default_bounds.penny.xy];
   end
end
for i=1:double(N_mctigue)
    if exist('mctigue_bounds','var') && length(mctigue_bounds.xy)>=5   bounds = [bounds; mctigue_bounds.xy(1:5,:)]; mctigue_bounds.xy(1:5,:)=[];  
       else                                                            bounds = [bounds; default_bounds.mctigue.xy];
   end
end
for i=1:double(N_yang)
    if exist('yang_bounds','var') && length(yang_bounds.xy) >= 8       bounds = [bounds; yang_bounds.xy(1:8,:)]; yang_bounds.xy(1:8,:)=[];     % FA 11/15  it was yang_bounds.xy>=5 which looks like a bug
       else                                                            bounds = [bounds; default_bounds.yang.xy];
   end
end
for i=1:double(N_squaredisloc)
    if exist('squaredisloc_bounds','var') && length(squaredisloc_bounds.xy) >= 9       bounds = [bounds; squaredisloc_bounds.xy(1:9,:)]; squaredisloc_bounds.xy(1:9,:)=[];  
       else                                                            bounds = [bounds; default_bounds.squaredisloc.xy];
   end
end
for i=1:double(N_multidisloc)
    multidisloc_len = 10 + length(multidislocopt.ind);
    if exist('multidisloc_bounds','var') && length(multidisloc_bounds.xy) >= multidisloc_len  bounds = [bounds; multidisloc_bounds.xy(1:multidisloc_len,:)]; multidisloc_bounds.xy(1:multidisloc_len,:)=[];  
       else                                                            bounds = [bounds; default_bounds.multidisloc.xy];
   end
end
for i=1:double(N_visco1d) 
    if exist('disloc_bounds','var') && length(disloc_bounds.xy) >= 10  bounds = [bounds; disloc_bounds.xy(1:10,:)]; disloc_bounds.xy(1:10,:)=[];  
       else                                                            bounds = [bounds; default_bounds.disloc.xy];
   end
end
for i=1:double(N_lockedandcreep)
    if exist('lockedandcreep_bounds','var') && length(lockedandcreep_bounds.xy) >= 5 bounds = [bounds; lockedandcreep_bounds.xy(1:5,:)]; lockedandcreep_bounds.xy(1:5,:)=[];  
       else                                                            bounds = [bounds; default_bounds.lockedandcreep.xy];
   end
end
for ni=1:double(N_peas)
    if exist('peas_bounds','var') && length(peas_bounds.xy) >=2 bounds = [bounds; peas_bounds.xy] ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  identify free and fixed parameters, modify bounds, parnames according to SARmul, etc  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linearind        = zeros(0,1);

freeind          = find( bounds(:,1) ~= bounds(:,2) );  
fixind           = find( bounds(:,1) == bounds(:,2) & bounds(:,1)~=999 );
linearind        = find( bounds(:,1) == 999         & bounds(:,2)==999 );
fixpar           = bounds(fixind,1);

remove           = [fixind;linearind];
bounds(remove,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% identify linear sources  (FA 7/2008: since all sources are linear in strength this should work by setting linearind=ones(1,N_sources)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linearsourceind  = zeros(1,N_sources);
tmplinearind     = linearind; 
i_source         = 1;
tmplinearind     = zeros(1,max([linearind;fixind;freeind]));
tmplinearind(linearind)=1;

for i=1:double(N_disloc) 
    ind=1:10; linearsourceind(i_source) = sum(tmplinearind(ind)); tmplinearind(ind)=[]; i_source=i_source+1;
end
for i=1:double(N_fault) 
    ind=1:9; linearsourceind(i_source) = sum(tmplinearind(ind)); tmplinearind(ind)=[]; i_source=i_source+1;
end

for i=1:double(N_mogi) 
    ind=1:4;  linearsourceind(i_source) = sum(tmplinearind(ind)); tmplinearind(ind)=[]; i_source=i_source+1;
end

for i=1:double(N_penny) 
    ind=1:5;  linearsourceind(i_source) = sum(tmplinearind(ind)); tmplinearind(ind)=[]; i_source=i_source+1;
end

for i=1:double(N_mctigue) 
    ind=1:5;  linearsourceind(i_source) = sum(tmplinearind(ind)); tmplinearind(ind)=[]; i_source=i_source+1;
end

for i=1:double(N_yang)
    ind=1:8;  linearsourceind(i_source) = sum(tmplinearind(ind)); tmplinearind(ind)=[]; i_source=i_source+1;
end

for i=1:double(N_multidisloc)
    ind=1:multidisloc_len;  linearsourceind(i_source) = sum(tmplinearind(ind)); tmplinearind(ind)=[]; i_source=i_source+1;
end

for i=1:double(N_lockedandcreep)
    ind=1:5;  linearsourceind(i_source) = sum(tmplinearind(ind)); tmplinearind(ind)=[]; i_source=i_source+1;
end

for i=1:double(N_peas)
    ind=1:6;  linearsourceind(i_source) = sum(tmplinearind(ind)); tmplinearind(ind)=[]; i_source=i_source+1;
end

if max(linearsourceind)>=2                                                                           
    errordlg('Only one linear parameter is allowed for each source. Sorry.'); error('-- exiting')  
end

%nonlinearsourceind = ~linearind;                                                          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  collect the inverse options 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pardim                          = max([max(fixind) max(freeind) max(linearind)]);     % FA 7/08 ParType for summary  

inverseopt.fixpar               = fixpar; 
inverseopt.fixind               = fixind; 
inverseopt.freeind              = freeind;
inverseopt.linearind            = linearind;
inverseopt.bounds               = bounds;

inverseopt.ParType              = cell(1,pardim);
inverseopt.ParType([fixind])    = {deal('fix')};
inverseopt.ParType([freeind])   = {deal('noli')};
inverseopt.ParType([linearind]) = {deal('lin')};

[inverseopt]                    = ModifyBoundsAndParNames(dataset,inverseopt,modelopt); 
inverseopt                      = orderfields(inverseopt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  collect the options that are actually used in the objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
objfuncopt.objfunc              = objfunc;

objfuncopt.modelopt             = modelopt;
objfuncopt.fixpar               = inverseopt.fixpar; 
objfuncopt.fixind               = inverseopt.fixind;
objfuncopt.freeind              = inverseopt.freeind;
objfuncopt.linearind            = inverseopt.linearind;
objfuncopt.linearsourceind      = logical(linearsourceind);

objfuncopt.PhaseRamp            = inverseopt.PhaseRamp;
objfuncopt.FactorLin            = inverseopt.FactorLin;
objfuncopt.FactorNonLin         = inverseopt.FactorNonLin;
objfuncopt.FactorNonLinDelta    = inverseopt.FactorNonLinDelta;

inverseopt.objfuncopt           = objfuncopt;
