function [multidislocopt]  = initialize_multidislocopt(multidislocopt,par_multidisloc)
%initialize_multidislocopt   - converts multidisloc option string into structure
 
if exist('par_multidisloc','var')
      if length(par_multidisloc)-10 ~= length(multidislocopt)-1
          errordlg('too much or too few entires in multidilocopt');error('-- exiting')
      end
  end
%
% first check whether argument is structure, contains an 'x' (e.g. 1x2) 
%
  if isstruct(multidislocopt) logmessage('multidislocopt is already a structure'); return; end
  if ~multidislocopt multidislocopt = '1x1'; end      % FA 2/2010   get some errors if this is not kept from version_1_0_3
%
% convert into cell array and extracts parameters from options cell array
%
  if ~iscell(multidislocopt) multidislocopt = str2list(multidislocopt); end

  all_names          = {multidislocopt{1:end}}; 
  [ind,parind]       = deal([]);
  [tmpind,tmpparind] = deal(zeros(size(all_names,2),1));
  
  for j=1:length(all_names)
      tmp               =   strmatch(all_names{j},{'Len' 'Wid' 'Dep' 'Dip' 'Strike' 'xE' 'yE' 'ss' 'ds' 'op' 'HorzOff' 'VertOff'});
      tmppar            =   strmatch(all_names{j},{'Len' 'Wid' 'Dep' 'Dip' 'Strike' 'xE' 'yE' 'ss' 'ds' 'op'});
      if ~isempty(tmp)      tmpind(j)    = tmp;    tmp    = []; end
      if ~isempty(tmppar)   
          tmpparind(j)  = tmppar; 
          tmppar        = []; 
      end
  end
  
  indVertOff = strmatch('VertOff',all_names);
  indHorzOff = strmatch('HorzOff',all_names);

  if length(indVertOff)~=length(indHorzOff) errordlg('multidislocopt not permissible; Need to specify HorzOff and VertOff for each dislocation'); error('--exiting'); end 

  N_disloc  = length(indVertOff) + 1;
  
  all_names = all_names(find(tmpind   ~=0));
  par_names = all_names(find(tmpparind~=0));
  
  tmpind   (find(tmpind   ==0))=[];
  tmpparind(find(tmpparind==0))=[];
  ind       = tmpind';
    
  for j=1:N_disloc-1;                        % FA 3/2010. used to have nc instead of N_disloc-1 [junk,nc] = mode(tmpparind);    % number of occurrences of most frequent parameter. 
       parind                             = [parind tmpparind(find(unique(tmpparind)))' + j*10];
       tmpparind(find(unique(tmpparind))) = [];
  end
   
  if ~isempty(par_names) && sum(parind(:))==0 errordlg([par_names{:} ' not allowed in multidislocopt. Use Len Wid Dep Dip Strike']); error('--exiting'); end

  new_parind                       = [1:10 parind];                   % add the index (1:10) from the original dislocation not included in all_names
  orig_parind                      = [1:10 10 + find(ind<=10)];
  orig_indVertOff                  = indVertOff + 10;
  orig_indHorzOff                  = indHorzOff + 10;
  
  clear multidislocopt;
  multidislocopt.all_names         = all_names;
  multidislocopt.ind               = ind;
  multidislocopt.orig_indVertOff   = orig_indVertOff;
  multidislocopt.orig_indHorzOff   = orig_indHorzOff;
  multidislocopt.orig_parind       = orig_parind;
  multidislocopt.new_parind        = new_parind;
  multidislocopt.N_disloc          = N_disloc;
  
  multidislocopt                   = orderfields(multidislocopt);
