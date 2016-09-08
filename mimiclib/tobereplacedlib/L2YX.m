function [Y,X]=L2YX(inputf,lenght)

% Gives 2D coordinates of the position in a row

X=ceil(inputf/lenght);
Y=inputf-((X-1)*lenght);