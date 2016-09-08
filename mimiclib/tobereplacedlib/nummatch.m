function [Indices]=nummatch(Nums1,Nums2)

% Find common value in two 1-column files, returns indices of value of vector 1, verctor 2 and common valuesvalue
% 
%               [Indices]=nummatch(Nums1,Nums2)
%
% N. Gourmelen
%

Indices=[];
p=0;
for i=1:size(Nums1,1)

 d=find(Nums2==Nums1(i));

 if d~=0
 
   p=p+1;
   Indices(p,:)=[i d Nums1(i)];

 end

end

if p==0

  p=1;
  Indices(p,:)=[-9999 -9999 -9999];

end

 

   
