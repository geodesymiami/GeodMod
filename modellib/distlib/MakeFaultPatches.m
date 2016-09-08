function [patches,nHorzEl,nVertEl,ExtendSize]=MakeFaultPatches(par,opt,direction)
%MakeFaultPatches             - Generates fault patches from dislocation
%
%  usage:  [modelopt]=MakeFaultPatches(distribopt,modelopt)
%
%  FA 7/2008 based on Sjonni Jonsson's code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  all default options are set in InverseDistrib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end
logmessage(sprintf('[]=%s(%s)',mfilename,inputname(1)));

if ~exist('direction','var') direction=1; end
switch direction
case{1}    
    
    if ~PatchSize  PatchSize = [par(1) par(2)]; end                   %set to length,width of original dislocation if 'off'

    tmp_disloc = par;

    if ExtendToSurface tmp_disloc      = Extend2Surface(tmp_disloc); end  % Extend Input dislocation to surface
    if ~ExtendSize     ExtendSize      = tmp_disloc(1:2);            end
    tmp_disloc(1:2) =  ExtendSize;
   
    nHorzEl  =  round( ExtendSize(1)/PatchSize(1) );
    nVertEl  =  round( ExtendSize(2)/PatchSize(2) );
    nEl      = nHorzEl * nVertEl;

    patches = patchfault(tmp_disloc(1:7)',nHorzEl,nVertEl);
    
case{-1}

otherwise
        errordlg('direction must be 1 or -1'); error('user error -- exiting');
end    
        

  
