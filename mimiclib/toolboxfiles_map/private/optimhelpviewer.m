function optimhelpviewer()
% OPTIMHELPVIEWER is a helper file for Optimtool. 

%   Copyright 2006-2007 The MathWorks, Inc. 
%   $Revision: 1.1.6.2 $

mapfilename = [docroot '/toolbox/optim/helptargets.map'];
try
    helpview(mapfilename, 'optimtool');
catch
    message = sprintf('Unable to display help for Optimtool');
    errordlg(message);
end
