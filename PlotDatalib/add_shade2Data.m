
function [igram]=add_shade2Data(Igram,ShadeFile,lopt)

%Prepare dem  and save it in the igram structure
%
%   Igram           = igram structure
%   ShadeFile           = DEM (.jpg only)
%

%   Noel Gourmelen October 2005
%   Falk Amelung 15 November
%                 Now allows lopt
%                 putting shade into igram(1).shade                               
%  FA Sep 2006, introduced use of fileparts funtion
%  FA Oct 2006, now wants the *jpg.rsc file consistent with the roi_pac names
%  FA Aug 2007, replaced ReadRscfile by ReadKeywordfile

% Read interfero

%  Igram=load(inigram);
  X_FIRSTi=Igram(1).x_first;Y_FIRSTi=Igram(1).y_first;
  X_STEPi=Igram(1).x_step;Y_STEPi=Igram(1).y_step;
  WIDTHi=size(Igram(1).data,2);FILE_LENGTHi=size(Igram(1).data,1);

  if isfield(Igram,'shade')
    if size(Igram(length(Igram)).shade,1) > 0
      
      error(sprintf('Exiting --- shade file exists'));
    end
  end

% Read dem 
 [pathstr, name, ext] = fileparts(ShadeFile);
 rscdem  = fullfile (pathstr , [name  ext '.rsc'])

  S=ReadKeywordfile(rscdem,' ');Sf=FilterKeyword(S,'dem');
  X_FIRSTd = S.X_FIRST;  Y_FIRSTd = S.Y_FIRST;
  X_STEPd  = S.X_STEP;    Y_STEPd  = S.Y_STEP;
  WIDTHd   = S.WIDTH;     FILE_LENGTHd=S.FILE_LENGTH;
 
  if length(isalpha(X_FIRSTd)) ~= 1  X_FIRSTd=str2num(X_FIRSTd) ; end
  if length(isalpha(Y_FIRSTd)) ~= 1  Y_FIRSTd=str2num(Y_FIRSTd) ; end
  if length(isalpha(X_STEPd)) ~= 1  X_STEPd=str2num(X_STEPd) ; end
  if length(isalpha(Y_STEPd)) ~= 1  Y_STEPd=str2num(Y_STEPd) ; end

  fid=fopen(ShadeFile,'r');
  [dem]=imread(ShadeFile,'jpeg') ;
  if  size(dem,3) == 1 
     tmp=zeros(size(dem,1),size(dem,2),3);
     tmp(:,:,1)=dem; 
     tmp(:,:,2)=dem; 
     tmp(:,:,3)=dem;
     dem=tmp ;
  end

  dem=double(dem);
  if exist('lopt','var')  && isfield(lopt,'subset')
      dems=Sf;dems.data=dem;
      dems=resize_igram(dems,lopt);
      dem=dems.data;
%      X_FIRSTd=dems;
%      Y_FIRSTd=Y_FIRSTd+lopt.subset.ji(1)*Y_STEPd;
  end

   dem(:,:,1)=dem(:,:,1)/256;dem(:,:,2)=dem(:,:,2)/256;dem(:,:,3)=dem(:,:,3)/256;

% Resize dem if diff pixel size

  factor=X_STEPd/X_STEPi;

  if factor < 0.99 | factor > 1.01
    factor
    dem=resizem(dem,factor);
    WIDTHd=size(dem,2);FILE_LENGTHd=size(dem,1);
  end


% Size the Dem to the Interfero

 % XoffF=(X_FIRSTi-X_FIRSTd)/X_STEPd+1;XoffE=XoffF+WIDTHi-1;
 % YoffF=(Y_FIRSTi-Y_FIRSTd)/Y_STEPd+1;YoffE=YoffF+FILE_LENGTHi-1;

 %DemRes=dem(round(YoffF):round(YoffE),round(XoffF):round(XoffE),:);

% Save the shade into igram

igram=Igram;
igram(length(igram)).shade=dem;
igram(length(igram)).x_first=igram(1).x_first;
igram(length(igram)).y_first=igram(1).y_first;
igram(length(igram)).x_step=igram(1).x_step;
igram(length(igram)).y_step=igram(1).y_step;

%save(inigram,'igram')

