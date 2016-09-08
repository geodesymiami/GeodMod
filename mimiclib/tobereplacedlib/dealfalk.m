function varargout = dealfalk(varargin)
%DEAL Deal inputs to outputs.
% same as deal but allows for variable output
% e.g.
% [D1,D2,D3,D4,D5]=deal(dataset.Ndata) gives an error
% if dataset.Ndata has only 4 array elements but 
% [D1,D2,D3,D4,D5]=dealfalk(dataset.Ndata) 
% returns D_5=[] or D_5=0 for logical quantities

% commented out some lines on Sep 12 05 because it did not work with only one input argument
%if nargin==1,
%  varargout = varargin(ones(1,max(1,nargout)));
%else
  tmp=cell(1,nargout);
  for i=1:nargin
      tmp{i} = varargin{i};
  end
  varargout = tmp;
  if islogical(varargout{1})          % assign false to empty fields
     for i=1:length(varargout)
         if isempty(varargout{i}) varargout{i}=false; end ;
     end
  end
%end
