function b=nanmedfilt2(a,m,n,maxnan)
% nanmedfilt2 - median filter ignoring NANs
%
% usage: newdata=nanmedfilt2(data,m,n,maxnan)
% 
% input:   data:    image to filter
%          m,n:     array size for median filtering
%          maxnan:  maximumm number of nan pixel allowed for each median computation
%                   (n*m fills wholes)                                       
%
% output:  newdata:  median filtered data
%
% FA 18 June 2001
%

if nargin==0, help nanmedfilt2; end; 
b=a;
mhalf=floor(m/2);
nhalf=floor(n/2);
for i=mhalf+1:size(a,1)-mhalf
for j=nhalf+1:size(a,2)-nhalf
    tmp=ones(1,n*m);
    tmp1 = a(i-mhalf:i+mhalf,j-nhalf:j+nhalf) ;
    tmp(:) = a(i-mhalf:i+mhalf,j-nhalf:j+nhalf) ;
	tmp(isnan(tmp))=[];
	if m*n-length(tmp) >= maxnan , tmp=nan; end          % check how many nan's
    b(i,j)=median(median( tmp ));
end
end

