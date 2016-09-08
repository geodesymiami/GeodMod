function [outvect] = getPolygon(data,opt)

%
%  Gather indices of an image by mouse clic.
%
%   [outvect] = getPolygon(data)
%
%
%  opt:
%
%     marKers :  [y y y;x x x] location of interest.
%
%
% N. Gourmelen, February 2008
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultopt=struct(   'marKers'                            ,        'NI'      ,             ...
                     'plotFigure'                         ,        'on'      ,             ...
                     'subsetArea'                         ,	       'off'     )              ;
    
if ~exist('opt','var')  [opt]=read_options_fromfile([mfilename '.min'],[]); end ;
[opt]=process_defaultoptions(opt,defaultopt);  
f=fieldnames(opt);  for i=1:length(f)  eval([char(f{i}) '= opt.(f{i}) ;' ]);  end

%%

ni = 1 ;  ny = 1 ;  testval = 1;

if plotFigure hh = figure ; set(hh,'Position',[0 100 600 450]) ;  imagesc(data) ;
    
    if length(marKers)~=1 ;  hold on ;  plot(marKers(2,:),marKers(1,:),'k^') ;  end
    
    while testval == 1
        
        user_in = input('Click yes when zoom, coloscale, ... is satisfactory [Y]: ', 's');
        if strcmp(lower(user_in),'y')  testval = 0 ; end
        
    end
    
end

testval = 1 ;    ['Select polygon:']

while testval == 1
    [li,Li,testvaln] = ginput(1) ;
    outvect(ny).poly(:,ni)   = round([Li;li]) ;
    hold on ; plot(outvect(ny).poly(2,:),outvect(ny).poly(1,:),'w-') ;
    ni = ni+1;

    if testvaln==3
        user_in = input('More polygon? Y/N [N]: ', 's');
        if isempty(user_in) user_in = 'N' ;  end
        if strcmp(lower(user_in),'n')  testval = 0 ;
        else testval = 1 ;  ny = ny +1 ;  ni = 1 ;  end
    end
end

