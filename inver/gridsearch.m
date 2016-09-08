function [mhat,misfit,list]=gridsearch(FUN,bounds,gridopt,varargin)
%gridsearch     [misfit,list]=gridsearch(FUN,bounds,gridopt,x1,x2...,xn)
%
%INPUTS:
%
%'FUN' specifies the objective function.  This function should accept an input model vector
%
%'bounds' specifies the upper and lower limits that each model parameter can take on.  This matrix
%must have as many rows as model parameters and two columns.
%
%'gridopt' specifies number of gridpoints for each parameter [if empty use 2]:
%
% TODO: (FA, noticed Feb 2008) there should be an option that gridopt is a vector with the grid spacing for each parameter
% TODO: (FA, noticed Feb 2008) list should be generated with a loop
%

%Check argument syntax

if nargin<2
	error('Usage: [mhat,F,model,energy,count]=gridseach(FUN,bounds,OPTIONS,varargin)')
end

if size(bounds,2)~=2
	error('Second argument must be an nx2 matrix of parameter bounds, where n is the number of parameters.');
end

%Check gridopt

	if isempty(gridopt)
		gridopt=2*ones(size(bounds,1),1);
	elseif length(gridopt)==1
		gridopt=gridopt*ones(size(bounds,1),1);
	end

%Check bounds to make sure they're ok

	if max(bounds(:,1)>bounds(:,2))
		error('All the values in the first column of bounds must be less than those in the second.');
	end

del=abs((bounds(:,1)-bounds(:,2)))./gridopt;
del(find(del==0))=1;
length(bounds(:,1));

switch length(bounds(:,1))

case{1}
  a1=bounds(1,1):del(1):bounds(1,2);
  [b1]=ndgrid(a1);
  list=[b1(:)];
case{2}                                      
  a1=bounds(1,1):del(1):bounds(1,2);
  a2=bounds(2,1):del(2):bounds(2,2);
  [b1,b2]=ndgrid(a1,a2);
  list=[b1(:),b2(:)];
case{3}                                      
  a1=bounds(1,1):del(1):bounds(1,2);
  a2=bounds(2,1):del(2):bounds(2,2);
  a3=bounds(3,1):del(3):bounds(3,2);
  [b1,b2,b3]=ndgrid(a1,a2,a3);
  list=[b1(:),b2(:),b3(:)];
case{3}                                      
  a1=bounds(1,1):del(1):bounds(1,2);
  a2=bounds(2,1):del(2):bounds(2,2);
  a3=bounds(3,1):del(3):bounds(3,2);
  [b1,b2,b3]=ndgrid(a1,a2,a3);
  list=[b1(:),b2(:),b3(:)];
case{4}                                      
  a1=bounds(1,1):del(1):bounds(1,2);
  a2=bounds(2,1):del(2):bounds(2,2);
  a3=bounds(3,1):del(3):bounds(3,2);
  a4=bounds(4,1):del(4):bounds(4,2);
  [b1,b2,b3,b4]=ndgrid(a1,a2,a3,a4);
  list=[b1(:),b2(:),b3(:),b4(:)];
case{5}                                      % Mogi
  a1=bounds(1,1):del(1):bounds(1,2);
  a2=bounds(2,1):del(2):bounds(2,2);
  a3=bounds(3,1):del(3):bounds(3,2);
  a4=bounds(4,1):del(4):bounds(4,2);
  a5=bounds(5,1):del(5):bounds(5,2);
  [b1,b2,b3,b4,b5]=ndgrid(a1,a2,a3,a4,a5);
  list=[b1(:),b2(:),b3(:),b4(:),b5(:)];
case{6}                                      % Mogi
  a1=bounds(1,1):del(1):bounds(1,2);
  a2=bounds(2,1):del(2):bounds(2,2);
  a3=bounds(3,1):del(3):bounds(3,2);
  a4=bounds(4,1):del(4):bounds(4,2);
  a5=bounds(5,1):del(5):bounds(5,2);
  a6=bounds(6,1):del(6):bounds(6,2);
  [b1,b2,b3,b4,b5,b6]=ndgrid(a1,a2,a3,a4,a5,a6);
  list=[b1(:),b2(:),b3(:),b4(:),b5(:),b6(:)];
case{7}                                      % Mogi
  a1=bounds(1,1):del(1):bounds(1,2);
  a2=bounds(2,1):del(2):bounds(2,2);
  a3=bounds(3,1):del(3):bounds(3,2);
  a4=bounds(4,1):del(4):bounds(4,2);
  a5=bounds(5,1):del(5):bounds(5,2);
  a6=bounds(6,1):del(6):bounds(6,2);
  a7=bounds(7,1):del(7):bounds(7,2);
  [b1,b2,b3,b4,b5,b6,b7]=ndgrid(a1,a2,a3,a4,a5,a6,a7);
  list=[b1(:),b2(:),b3(:),b4(:),b5(:),b6(:),b7(:)];
case{8}                                    % Yang
  a1=bounds(1,1):del(1):bounds(1,2);
  a2=bounds(2,1):del(2):bounds(2,2);
  a3=bounds(3,1):del(3):bounds(3,2);
  a4=bounds(4,1):del(4):bounds(4,2);
  a5=bounds(5,1):del(5):bounds(5,2);
  a6=bounds(6,1):del(6):bounds(6,2);
  a7=bounds(7,1):del(7):bounds(7,2);
  a8=bounds(8,1):del(8):bounds(8,2);
  [b1,b2,b3,b4,b5,b6,b7,b8]=ndgrid(a1,a2,a3,a4,a5,a6,a7,a8);
  list=[b1(:),b2(:),b3(:),b4(:),b5(:),b6(:),b7(:),b8(:)];
case{9}                                    % Yang
  a1=bounds(1,1):del(1):bounds(1,2);
  a2=bounds(2,1):del(2):bounds(2,2);
  a3=bounds(3,1):del(3):bounds(3,2);
  a4=bounds(4,1):del(4):bounds(4,2);
  a5=bounds(5,1):del(5):bounds(5,2);
  a6=bounds(6,1):del(6):bounds(6,2);
  a7=bounds(7,1):del(7):bounds(7,2);
  a8=bounds(8,1):del(8):bounds(8,2);
  a9=bounds(9,1):del(9):bounds(9,2);
  [b1,b2,b3,b4,b5,b6,b7,b8,b9]=ndgrid(a1,a2,a3,a4,a5,a6,a7,a8,a9);
  list=[b1(:),b2(:),b3(:),b4(:),b5(:),b6(:),b7(:),b8(:),b9(:)];
case{10}                                    % Yang
  a1=bounds(1,1):del(1):bounds(1,2);
  a2=bounds(2,1):del(2):bounds(2,2);
  a3=bounds(3,1):del(3):bounds(3,2);
  a4=bounds(4,1):del(4):bounds(4,2);
  a5=bounds(5,1):del(5):bounds(5,2);
  a6=bounds(6,1):del(6):bounds(6,2);
  a7=bounds(7,1):del(7):bounds(7,2);
  a8=bounds(8,1):del(8):bounds(8,2);
  a9=bounds(9,1):del(9):bounds(9,2);
  a10=bounds(10,1):del(10):bounds(10,2);
  [b1,b2,b3,b4,b5,b6,b7,b8,b9,b10]=ndgrid(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10);
  list=[b1(:),b2(:),b3(:),b4(:),b5(:),b6(:),b7(:),b8(:),b9(:),b10(:)];
case{11}                                    % Yang
  a1=bounds(1,1):del(1):bounds(1,2);
  a2=bounds(2,1):del(2):bounds(2,2);
  a3=bounds(3,1):del(3):bounds(3,2);
  a4=bounds(4,1):del(4):bounds(4,2);
  a5=bounds(5,1):del(5):bounds(5,2);
  a6=bounds(6,1):del(6):bounds(6,2);
  a7=bounds(7,1):del(7):bounds(7,2);
  a8=bounds(8,1):del(8):bounds(8,2);
  a9=bounds(9,1):del(9):bounds(9,2);
  a10=bounds(10,1):del(10):bounds(10,2);
  a11=bounds(11,1):del(11):bounds(11,2);
  [b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11]=ndgrid(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11);
  list=[b1(:),b2(:),b3(:),b4(:),b5(:),b6(:),b7(:),b8(:),b9(:),b10(:),b11(:)];

case{12}
  a1=bounds(1,1):del(1):bounds(1,2);
  a2=bounds(2,1):del(2):bounds(2,2);
  a3=bounds(3,1):del(3):bounds(3,2);
  a4=bounds(4,1):del(4):bounds(4,2);
  a5=bounds(5,1):del(5):bounds(5,2);
  a6=bounds(6,1):del(6):bounds(6,2);
  a7=bounds(7,1):del(7):bounds(7,2);
  a8=bounds(8,1):del(8):bounds(8,2);
  a9=bounds(9,1):del(9):bounds(9,2);
  a10=bounds(10,1):del(10):bounds(10,2);
  a11=bounds(11,1):del(11):bounds(11,2);
  a12=bounds(12,1):del(12):bounds(12,2);
  [b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12]=ndgrid(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12);
  list=[b1(:),b2(:),b3(:),b4(:),b5(:),b6(:),b7(:),b8(:),b9(:),b10(:),b11(:),b12(:)];
case{13}
  a1=bounds(1,1):del(1):bounds(1,2);
  a2=bounds(2,1):del(2):bounds(2,2);
  a3=bounds(3,1):del(3):bounds(3,2);
  a4=bounds(4,1):del(4):bounds(4,2);
  a5=bounds(5,1):del(5):bounds(5,2);
  a6=bounds(6,1):del(6):bounds(6,2);
  a7=bounds(7,1):del(7):bounds(7,2);
  a8=bounds(8,1):del(8):bounds(8,2);
  a9=bounds(9,1):del(9):bounds(9,2);
  a10=bounds(10,1):del(10):bounds(10,2);
  a11=bounds(11,1):del(11):bounds(11,2);
  a12=bounds(12,1):del(12):bounds(12,2);
  a13=bounds(13,1):del(13):bounds(13,2);
  [b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13]=ndgrid(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13);
  list=[b1(:),b2(:),b3(:),b4(:),b5(:),b6(:),b7(:),b8(:),b9(:),b10(:),b11(:),b12(:),b13(:)];
end

list  =list';
misfit=zeros(size(list,1),1);

tic
for i=1:size(list,2)
	f=feval(FUN,list(:,i),varargin{:});
	misfit(i)=f;
end

mhat=list(:,find(misfit==min(misfit)));

tt=toc;
sprintf('%d function calls: %f seconds or %f hrs', size(list,1),  tt, tt/3600)


