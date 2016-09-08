function [outvect] = get_polygon(data)

%
%  Gather indices of an image by mouse clic.
%
%   [outvect] = get_polygon(data)
%
%
% N. Gourmelen, February 2008
%

ni = 1 ;  ny = 1 ;  testval = 1;

hh = figure ; set(hh,'Position',[0 100 600 450]) ;  imagesc(data) ;

while testval == 1
    [li,Li,testvaln] = ginput(1) ;
    outvect(ny).poly(:,ni)   = round([Li;li]) ;
    hold on ; plot(outvect(ny).poly(2,:),outvect(ny).poly(1,:),'w-') ;
    ni = ni+1;

    if testvaln==3
        user_in = input('More polygon? Y/N [N]: ', 's');
        if isempty(user_in) user_in = 'N' ;  end
        if strcmp(user_in,'N')  testval = 0 ;
        else testval = 1 ;  ny = ny +1 ;  ni = 1 ;  end
    end
end

