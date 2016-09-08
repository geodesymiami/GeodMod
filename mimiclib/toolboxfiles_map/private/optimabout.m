function optimabout()
%OPTIMABOUT helper that displays the About Box  

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2007/06/14 05:19:11 $

a = ver('optim');
str = sprintf(['Optimization Toolbox ' a.Version '\n',...
               'Copyright 1990-' a.Date(end-3:end) ' The MathWorks, Inc.']);
msgbox(str,'About Optimization Toolbox', 'modal');
