function [L]=YX2L(Y,X,File_length)

% Transform Y X indices to vector
%
% [L]=YX2L(Y,X,File_length)
%

L=(X-1)*File_length+Y;

