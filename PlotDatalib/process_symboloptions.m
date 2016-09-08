function [symbols]=process_symboloptions(plotopt,symbolopts);
%process_symbolsoptions     - generates symbol structure ; 
%
%usage: [symbols]=process_symboloptions(symbols,default_symbols);
%
%Input:   plotopt            structure containing fields for which symbols are generated
%         symbolopts         structure containing desired symbols          
%
%Output:  symbols    new structure with symbols  
%
% NOTE: I wanted to use a symbopt structure of the form
%       symbopt.Faults.LineWidth  = 5 ;
%       symbopt.Triangles.MarkerSize = 10 ;
% which can be easily converted into the symbols structure but I don't know how to address the default value ('w' for symbols.Faults)
% because there is no corresponding FieldName as for example LineWidth. Maybe there is a name that the plot function understands ?
% If not we could call it symbopt.Faults.Default='w' and have process_symbopt treating the Default field differently
%FA, Nov 2005  modified from process_options.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Default plot symbols %%%%%%%%%
default_symbols.Faults       =  {'w'}    ;   
default_symbols.Roads        =  {'w'}    ;   
default_symbols.Vectors      =  {'k'}    ;   
default_symbols.Triangles    =  {'k^'}   ;   
default_symbols.VectorsBlack =  {'k'}    ;   
default_symbols.VectorsRed   =  {'r'}    ;   
        
default_symbols.Quakes.mag_ranges   = [2 3 9] ;            
default_symbols.Quakes.mag_symbols  = { {'MarkerSize',10} {'MarkerSize',20} } ;    % one cell less than ranges
default_symbols.Quakes.dep_ranges = [0 8 15 100] ;       
default_symbols.Quakes.dep_symbols= {'ro' 'bo' 'yo'} ;                         % one cell less than ranges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%  generate empty symbols structure (e.g.: symbols.Faults={} )
%

f=fieldnames(plotopt) ;
for i=1:length(f)
    eval(['symbols.(f{i}) = {} ;' ]) ;
end

%
%   assign default symbols  if given                       %
%

f=fieldnames(default_symbols) ;  
for i=1:length(f)
    symbols.(f{i})=default_symbols.(f{i}) ; 
end

%
%   assign symbopts if given                       %
%

if ~isempty(symbolopts) && isstruct(symbolopts)
   f=fieldnames(symbolopts) ;  
   for i=1:length(f)
       if isstruct(symbols.(f{i}))
          g=fieldnames(symbols.(f{i})) ;  
          for ii=1:length(g)
              if  isfield(symbolopts.(f{i}),g{ii}) && ~isempty(symbolopts.(f{i}).(g{ii}))
                  symbols.(f{i}).(g{ii})=symbolopts.(f{i}).(g{ii});
              end
          end
       else
          symbols.(f{i})=symbolopts.(f{i}) ;
       end
   end
end

%
%  display options on screen
%

f=fieldnames(symbols) ;  
for i=1:length(f)
    if  isfield(symbols,(f{i})) && ~isempty(symbols.(f{i}))
        str1=['symbols.' f{i} ':']; str2=[symbols.(f{i})];
        %disp(str1),disp(str2);
    end
end
