function [block] = BlockDiag(A1,A2,A3,A4,A5,A6,A7,A8)
%   BlockDiag      - Creates a block diagonal matrix
%
% Input:        A1, A2, A3.... up to 8 matrices
%
% Output:       block = a block diagona matrix on the form
%
%                               | A1  0  0 ... |
%                       block = |  0 A2  0 ... |
%                               |  0  0 A3 ... |

block = A1; 
  
  if nargin>1; [lin,col]=size(block); [li2,co2]=size(A2); block = [block zeros(lin,co2); zeros(li2,col) A2]; end
  if nargin>2; [lin,col]=size(block); [li2,co2]=size(A3); block = [block zeros(lin,co2); zeros(li2,col) A3]; end
  if nargin>3; [lin,col]=size(block); [li2,co2]=size(A4); block = [block zeros(lin,co2); zeros(li2,col) A4]; end
  if nargin>4; [lin,col]=size(block); [li2,co2]=size(A5); block = [block zeros(lin,co2); zeros(li2,col) A5]; end
  if nargin>5; [lin,col]=size(block); [li2,co2]=size(A6); block = [block zeros(lin,co2); zeros(li2,col) A6]; end
  if nargin>6; [lin,col]=size(block); [li2,co2]=size(A7); block = [block zeros(lin,co2); zeros(li2,col) A7]; end
  if nargin>7; [lin,col]=size(block); [li2,co2]=size(A8); block = [block zeros(lin,co2); zeros(li2,col) A8]; end
  if nargin>8; disp('Too many arguments'); end
 
 