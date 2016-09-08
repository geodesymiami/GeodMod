function [ind] = find_coord_index(modelopt);
%   finds_coord_index - find index vector for x,y coordinates of model parameters (for lola2xy conversion)
% usage:  [ind]=find_coord_index(objfuncopt);
%
% FA, May 2007 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind=[];
m_length=length(modelopt.ParNames);

N_disloc     = modelopt.N_disloc;      
N_fault      = modelopt.N_fault;      
N_mogi       = modelopt.N_mogi;        
N_penny      = modelopt.N_penny;       
N_mctigue    = modelopt.N_mctigue;     
N_yang       = modelopt.N_yang;        
N_multidisloc= modelopt.N_multidisloc;        
N_visco1d    = modelopt.N_visco1d;     
%N_lockedandcreep = modelopt.N_lockedandcreep   % FA Feb 2008: not needed because lockedandcreep does not have coordinates as modelparameters

if isfield(modelopt,'multidislocopt') multidisloc_len = 10+length(modelopt.multidislocopt.ind); else multidisloc_len=[]; end

%
% calculate the index of the first parameter of the first source of each type.
% The index depends on hown many sources are at lower indices.
% For example, for 2 Dislocations (each 10 parameters), 1 Mogi (4 par) and 2 Penny sources (each 5 par)
% the index of the first mctigue source is 35= 2*10 +1*4 + 2*5
% 
ifirst_disloc      = 1;
ifirst_fault       = 1;              
ifirst_mogi        = ifirst_disloc     +  N_disloc*10;
ifirst_penny       = ifirst_mogi       +  N_mogi*4;
ifirst_mctigue     = ifirst_penny      +  N_penny*5;
ifirst_yang        = ifirst_mctigue    +  N_mctigue*5;
ifirst_multidisloc = ifirst_yang       +  N_yang*8;
ifirst_visco1d     = ifirst_multidisloc +  multidisloc_len;

%
% now calculate the indices of the x,y or long,lat parameters ([6 7 16 17] for 2 dislocations)
%
if N_disloc >= 1      ind = [ind  ifirst_disloc+5  ifirst_disloc+6 ] ;  end
if N_disloc >= 2      ind = [ind  ifirst_disloc+15 ifirst_disloc+16] ;  end
if N_disloc >= 3      ind = [ind  ifirst_disloc+25 ifirst_disloc+26] ;  end
if N_disloc >= 4      ind = [ind  ifirst_disloc+35 ifirst_disloc+36] ;  end
if N_disloc >= 5      for i=5:N_disloc  ind = [ind  ifirst_disloc+(i-1)*10+5 ifirst_disloc+(i-1)*10+6]; end; end

if N_fault  >= 1      ind = [ind  ifirst_fault  ifirst_fault+1 ] ;  end

if N_mogi   >= 1 ind = [ind  ifirst_mogi      ifirst_mogi+1 ] ;  end
if N_mogi   >= 2 ind = [ind  ifirst_mogi+4    ifirst_mogi+5 ] ;  end
if N_mogi   >= 3 ind = [ind  ifirst_mogi+8    ifirst_mogi+9 ] ;  end

if N_penny       == 1 ind = [ind  ifirst_penny         ifirst_penny+1      ];  end
if N_mctigue     == 1 ind = [ind  ifirst_mctigue       ifirst_mctigue+1    ];  end
if N_yang        == 1 ind = [ind  ifirst_yang          ifirst_yang+1       ];  end
if N_multidisloc == 1 ind = [ind  ifirst_multidisloc+5 ifirst_multidisloc+6];  end
if N_visco1d     >= 1 ind = [ind  ifirst_visco1d+5     ifirst_visco1d+6    ];  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
