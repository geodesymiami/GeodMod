function [rangechange]=displacement2rangechange(u,radarlook)
%displacement2rangechange  - converts displacement into range change 
%
% INPUT:
%		u	             - displacement vector
%       radarlook        - unit vector in LOS direction            
%       
% OUTPUT:
%       rangechange
%
%

    udim  =size(u,2);

                      % the following works fine for udim=1
                      %rangechange = zeros(length(u)/3,udim);
                      %utmp        = u.*radarlook(:); 
                      %utmp3       = reshape(utmp,3,[]);
                      %rangechange = sum(utmp3)';
    
    rangechange    = zeros(length(u)/3,udim);

    utmp           = u.*repmat(radarlook(:),1,udim);      % multiply element-wise    (FA 7/2008)
    utmp3          = reshape(utmp,3,[]);                  % reshape into a 3xNdata arrary
    rangechangetmp = sum(utmp3)';                         % summ over the 3 components (scalar product)
    rangechange    = reshape(rangechangetmp,[],udim);     % reshape rangechange into one colum for each source
