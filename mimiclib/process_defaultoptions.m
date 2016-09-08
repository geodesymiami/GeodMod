function [optionsopt]=process_defaultoptions(optionsopt,defaultopt);
%process_defaultoptions     - assigns default options to empty fields in option structure
%
%usage: [optionsopt]=process_options(optionsopt,defaultopt);
%
%Input:   optionsopt    structure containing options (e.g. objfuncopt, can be empty)             
%         defaultopt    structure containing default options (can be a structure)                             
%         all fields 'off','on' are set to false,true
%
%Output:  optionsopt    new structure with options  
%
%Currently supports 1 nested structure (e.g. plotdataopt.colorbar.Location='off') 
%TODO:  need to modify to support 2 nested structures 
% (e.g. rates2threedfieldopt.plotdataopt.colorbaropt.Location= 'OutsideLowerRight') 
%
%FA, June 4, 2005  modified from process_options.m
%FA, Nov 6, 2006   now allows a field to contain a structure
%FA, May 2, 2007   now also sets 'off' and 'on' fields in substructures to false,true

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   set inverseopt                       %
f=fieldnames(defaultopt) ;  
if ~exist('optionsopt','var')  optionsopt=[]  ; end
for i=1:length(f)
    if isstruct(defaultopt.(f{i}))
       if ~isfield(optionsopt,f{i}) 
           optionsopt.(f{i})=defaultopt.(f{i}) ; 
       else
           g=fieldnames(defaultopt.(f{i})) ;  
           for j=1:length(g)
               if ~isfield(optionsopt.(f{i}),g{j}) 
                   optionsopt.(f{i}).(g{j})=defaultopt.(f{i}).(g{j}) ; 
               end
           end
       end
       optionsopt.(f{i})=orderfields(optionsopt.(f{i}));
    else
       if ~isfield(optionsopt,f{i}) optionsopt.(f{i})=defaultopt.(f{i}) ; end
    end
end

optionsopt=orderfields(optionsopt);

% now replace all 'off' by false and all 'on' by true
f=fieldnames(optionsopt) ;  
for i=1:length(f)
    if  strcmp('off',optionsopt.(f{i})) optionsopt.(f{i})=false ; end
    if  strcmp('on' ,optionsopt.(f{i})) optionsopt.(f{i})=true  ; end

    if isstruct(optionsopt.(f{i})) 
       g=fieldnames(optionsopt.(f{i})) ;  
       for j=1:length(g)
           if  strcmp('off',optionsopt.(f{i}).(g{j})) optionsopt.(f{i}).(g{j})=false  ; end
           if  strcmp('on' ,optionsopt.(f{i}).(g{j})) optionsopt.(f{i}).(g{j})=true   ; end
       end
    end
end
   

%% Code before May 2 2007 changes
%% now replace all 'off' by false
%f=fieldnames(optionsopt) ;  
%for i=1:length(f)
%    if ~isstruct(optionsopt.(f{i})) &  strcmp('off',optionsopt.(f{i})) optionsopt.(f{i})=false ; end
%end
%% now replace all 'on' by true
%f=fieldnames(optionsopt) ;  
%for i=1:length(f)
%    if ~isstruct(optionsopt.(f{i})) &  strcmp('on',optionsopt.(f{i})) optionsopt.(f{i})=true ; end
%end
