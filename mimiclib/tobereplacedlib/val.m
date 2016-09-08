function [result]=val(file,Lindices,Cindices)

%find values of matrix when indices are vectors

lim=length(Lindices);
result=ones(lim,1);

if size(Lindices,2)==2
    Cindices = Lindices(:,2) ;
    Lindices = Lindices(:,1) ;
end

for i=1:lim

 if Lindices(i) > size(file,1) | Cindices(i) > size(file,2)
  
  result(i,1)=NaN;

 else

  result(i,1)=file(Lindices(i),Cindices(i));
 
 end

end
