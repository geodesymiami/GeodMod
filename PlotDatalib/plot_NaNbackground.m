
function plot_NaNbackground(data,opt)

%
% Plot data containg NaNs. Only data will be used to scale the colorscale
%
%   data  :  data!
%   
%   opt:
%   
%       cmap             : Colormap 
%       background_color : Example -> [0] for white (default)
%
%
%   Noel Gourmelen - March 2009
%

defaultopt=struct(                                                       	      ...
    'cmap'           		                    ,        'off'      ,             ...
    'background_color'                              ,        'off'      )              ;
    

if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],[]); end ;
[opt] = process_defaultoptions(opt,defaultopt);  display(opt)
f=fieldnames(opt) ;
for i=1:length(f)
    eval([char(f{i}) '= opt.(f{i}) ;' ]) ;
end

z = data ;


NaNi = find(isnan(data)==1) ;
if ~cmap  cmap = colormap ;  end 
if ~background_color  background_color = [0] ;  end

cmap=[background_color;cmap] ;

if isempty(NaNi)  
	sprintf('Contains no NaNs') ;
else  
	data(NaNi(1)) = min(data(:)) - ( max(data(:)) -min(data(:)) ) / length(cmap) ;
end

imagesc(data) ; colormap(cmap)
