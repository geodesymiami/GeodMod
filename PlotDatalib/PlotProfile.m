function PlotProfile(Xdata,Ydata,opt)

%
% Plot Profiles
%
%  PlotProfile(Xdata,Ydata,opt)
%
%   inputs:
%
%       Xdata: X axis data vector (1*n)
%       
%       Ydata: Y axis data vector (can plot multiple vectors -> x*n)
%
%       opt can be:
%
%           .y_unit      
%
%           .profilecolor -> string (1*x) of color codes. Default is ['k';'r';'b';'g';'y';'m';'c']
%
%           .profilestyle -> string (1*x) of style codes. Default is ['.']
%
%
% N. Gourmelen. Feb. 2008
%

defaultopt=struct(                                                                          ...
        'DoIt'               ,                            'on '                ,            ...
        'profilecolor'       ,     ['k';'r';'b';'g';'y';'m';'c']               ,            ...
        'profilestyle'       ,                              '.'                ,            ...
        'profiletitle'       ,                   'Profile Plot'                ,            ...
        'y_unit'             ,                          'm'                    ,            ...
        'profileXlabel'      ,                         'Distance (km)'         ,            ...
        'profileYlabel'      ,                         'Displacement'          )             ;
    
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],defaultopt); end ;
[opt]=process_defaultoptions(opt,defaultopt);           %display(opt)
f=fieldnames(opt) ; for i=1:length(f) eval([char(f{i}) '= opt.(f{i}) ;' ]) ; end
if  ~DoIt                         return; end

%% Check input data

if size(Xdata,2) ~= size(Ydata,2)        ;  error('Error in PlotProfile, dimension of Xdata and Ydata should agree.')                 ;  end
if size(Ydata,1)  > size(profilecolor,1) ;  error('Error in PlotProfile, Number of specified color should equal number of profiles.') ;  end
    

%% Plot

figure;

for ni=1:size(Ydata,1)
    plotstyle = [profilecolor(ni),profilestyle];
    plot(Xdata,Ydata(ni,:),plotstyle);  
    hold on;
end

title(profiletitle) ;  xlabel(profileXlabel) ;  ylabel([profileYlabel,' (',y_unit,')']) ;





