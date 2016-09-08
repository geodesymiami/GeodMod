%Max value of array p1:
function [Maxf,minf]= maxmin (file); 

maxp1=max(file); %by columns
Maxf=max(maxp1);  

minp1=min(file); %by columns
minf=min(minp1);