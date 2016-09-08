function [modelopt,default_bounds]=InitializeModelopt(modelopt,basemap)
%InitializeModelopt  -  default model options and inversion bounds  
%
%usage:  [modelopt]                = InitializeModelopt(modelopt)
%        [modelopt]                = InitializeModelopt(modelopt,basemap)
%        [modelopt,default_bounds] = InitializeModelopt(modelopt,basemap)
%
%  Sets default model options. In the second case par_lola is calculated from par_xy, if given
%  ans vice-versa. In the third case the default bounds are calculated.
%  
%  Uses N_disloc=1 if no sources given.
%
%  The convention of source order is:
%      1. disloc
%      2. mogi
%      3. penny
%      4. mctigue
%      5. yang
%      6. multidisloc
%      7. visco1d
%      8. lockedandcreep
%
% TODO: At the end of the function need to clean out (rmfield) visco1dopt, etc 
% V1.0  Falk Amelung, June 2007, modified from MakeParnames.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process default options, and set variables to options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultopt=struct(                                                         ...
        'N_disloc'           ,        'off'                   ,            ...
        'N_fault'            ,        'off'                   ,            ...
        'N_mogi'             ,        'off'                   ,            ...
        'N_penny'            ,        'off'                   ,            ...
        'N_mctigue'          ,        'off'                   ,            ...
        'N_yang'             ,        'off'                   ,            ...
        'N_squaredisloc'     ,        'off'                   ,            ...
        'N_multidisloc'      ,        'off'                   ,            ...
        'N_visco1d'          ,        'off'                   ,            ...
        'N_lockedandcreep'   ,        'off'                   ,            ...
        'N_peas'             ,        'off'                   ,            ... 
        'N_sources'          ,        'off'                   ,            ...
        'Topo'               ,        'off'                   ,            ...
        'Layers'             ,        'off'                   ,            ...
        'multidislocopt'     ,        'off'                   ,            ...
        'Unit'               ,        'm'        )            ;
defaultopt.dislocopt=struct(                         ...                                                                                   
        'nu'                 ,         0.25      )            ;
%defaultopt.multidislocopt=struct(                    ...
%       'len'                ,        'off'                   ,           ...
%       'matrix'             ,        [1 1]                   ,           ...
%       'nu'                 ,         0.25      )            ;
defaultopt.visco1dopt=struct(                        ...                                                                                   
        'nu'                 ,         0.25      )            ;

[modelopt]                = process_defaultoptions(modelopt,defaultopt);  
[modelopt.multidislocopt] = initialize_multidislocopt(modelopt.multidislocopt);
f=fieldnames(modelopt) ; for i=1:length(f) eval([char(f{i}) '= modelopt.(f{i}) ;' ]) ; end
logmessage(mfilename);

N_sources = N_disloc+N_fault+N_mogi+N_penny+N_mctigue+N_yang+N_multidisloc+N_visco1d+N_lockedandcreep+N_peas+N_squaredisloc;

if N_sources==0  N_disloc=1 ; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% calculate bounds of area for default bounds from basemap if given %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   x_min   = 0 ;
   y_min   = 0 ;
   if exist('basemap','var')
      x_max  = basemap.x_posting*size(basemap.data,2);
      y_max  = basemap.y_posting*size(basemap.data,1);
   else
      x_max  = 50 ;
      y_max  = 50 ;
   end
   x_start = x_min + (x_max-x_min)/4 ;
   x_end   = x_max - (x_max-x_min)/4 ;
   y_start = y_min + (y_max-y_min)/4 ;
   y_end   = y_max - (y_max-y_min)/4 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% set format strings and default bounds %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disloc_names         = {'Len'  'Wid'    'Dep'   'Dip'   'Strike' 'xE'    'xN'    'ss'    'ds'    'op'   };   
disloc_namform       = {'%5s'  '%5s'    '%5s'   '%7s'   '%7s'    '%6s'   '%6s'   '%5s'   '%5s'   '%5s'  };   
disloc_varform       = {'%5.2f' '%5.2f' '%5.2f' '%7.2f' '%7.2f'  '%6.2f' '%6.2f' '%5.2f' '%5.2f' '%5.2f'};
disloc_bounds.xy(:,1)= [  2       2       2     -89.9     0.1    x_start y_start  0       0       999   ];
disloc_bounds.xy(:,2)= [ 10      10      10      89.9   179.9    x_end   y_end    0       0       999   ];

fault_names         = {'xE'      'xN'     'Strike' 'Dip'   'Rake'   'Len'   'Top'   'Bot'   'Slip'   };   
fault_namform       = {'%5s'     '%5s'    '%7s'    '%7s'   '%7s'    '%5s'   '%5s'   '%5s'   '%5s'    };   
fault_varform       = {'%5.2f'   '%5.2f'  '%7.2f'  '%7.2f' '%7.2f'  '%5.2f' '%5.2f' '%5.2f' '%5.2f'  };
fault_bounds.xy(:,1)= [  x_start y_start    0.1    -89.9    -179.9     2       2       5      0 ];
fault_bounds.xy(:,2)= [  x_end   y_end      179.9   89.9     179.9    10       5      10      0 ];

mogi_names          = {'xE'    'xN'    'Dep'   'Stren'};                                        
mogi_namform        = {'%5s'   '%5s'   '%5s'   '%5s'  };                                        
mogi_varform        = {'%5.2f' '%5.2f' '%5.2f' '%5.2f'};
mogi_bounds.xy(:,1) = [x_start y_start  2        999  ];
mogi_bounds.xy(:,2) = [x_end   y_end    10       999  ];

penny_names         = {'xE'    'xN'    'Dep'   'Rad'   'Stren'};                                        
penny_namform       = {'%5s'   '%5s'   '%5s'   '%5s'   '%5s'  };                                        
penny_varform       = {'%5.2f' '%5.2f' '%5.2f' '%5.3f' '%5.2f'};
penny_bounds.xy(:,1)= [x_start y_start  2        1       999  ];
penny_bounds.xy(:,2)= [x_end   y_end    10       2       999  ];

mctigue_names         = {'xE'    'xN'    'Dep'   'Rad'   'en'};                                        
mctigue_namform       = {'%5s'   '%5s'   '%5s'   '%5s'   '%5s'  }  ;                                        
mctigue_varform       = {'%5.2f' '%5.2f' '%5.2f' '%5.3f' '%5.2f'};
mctigue_bounds.xy(:,1)= [x_start y_start  2        1       999  ];
mctigue_bounds.xy(:,2)= [x_end   y_end    10       2       999  ];

yang_names          = {'xE'    'xN'    'Dep'   'Press' 'majAx' 'AxRatio' 'Strike' 'Plunge' };                                        
yang_namform        = {'%5s'   '%5s'   '%5s'   '%5s'   '%5s'   '%7s'     '%7s'    '%7s'    };                                        
yang_varform        = {'%5.2f' '%5.2f' '%5.2f' '%5.3f' '%5.2f' '%5.2f'   '%7.2f'  '%6.2f'  };
yang_bounds.xy(:,1) = [x_start y_start  2        999     1.0     0.5       0.1       0.01  ];
yang_bounds.xy(:,2) = [x_end   y_end    10       999     1.0     0.5     179.9       0.01  ];

squaredisloc_names         = {'Len'   'Dep'   'Dip'   'Strike' 'xE'    'xN'    'ss'    'ds'    'op'   };   
squaredisloc_namform       = {'%5s'   '%5s'   '%7s'   '%7s'    '%6s'   '%6s'   '%5s'   '%5s'   '%5s'  };   
squaredisloc_varform       = {'%5.2f' '%5.2f' '%7.2f' '%7.2f'  '%6.2f' '%6.2f' '%5.2f' '%5.2f' '%5.2f'};
squaredisloc_bounds.xy(:,1)= [  2       2     -89.9     0.1    x_start y_start  0       0       999   ];
squaredisloc_bounds.xy(:,2)= [ 10      10      89.9   179.9    x_end   y_end    0       0       999   ];

tmp_multidisloc_names         = {disloc_names{:}  'HorzOff' 'VertOff'};
tmp_multidisloc_namform       = {disloc_namform{:}  '%7s'   '%7s'    };
tmp_multidisloc_varform       = {disloc_varform{:}  '%7.2f' '%7.2f'  };
tmp_multidisloc_bounds.xy(:,1)= [disloc_bounds.xy(:,1)'  0     0     ]';
tmp_multidisloc_bounds.xy(:,2)= [disloc_bounds.xy(:,2)'  5     5     ]';

multidisloc_names             = {disloc_names{:}   tmp_multidisloc_names{multidislocopt.ind}};
multidisloc_namform           = {disloc_namform{:} tmp_multidisloc_namform{multidislocopt.ind}};
multidisloc_varform           = {disloc_varform{:} tmp_multidisloc_varform{multidislocopt.ind}};
multidisloc_bounds.xy         = [disloc_bounds.xy; tmp_multidisloc_bounds.xy(multidislocopt.ind,:)];

lockedandcreep_names          = {'ld'    'dsd'   'usd'  'loff'  'ffv'   };                                        
lockedandcreep_namform        = {'%5s'   '%5s'   '%5s'   '%5s'   '%8s'  };                                        
lockedandcreep_varform        = {'%5.2f' '%5.2f' '%5.2f' '%5.2f' '%8.5f'};
lockedandcreep_bounds.xy(:,1) = [   0      0       0      0       0     ];
lockedandcreep_bounds.xy(:,2) = [  10     10      10      0       0.005 ];

peas_names          = {'height'  'porePress'  'depth'  'poros' 'injRate' 'wellLength' 'mu' 'permeab'} ;   
peas_namform        = {'%6s'   '%10s'   '%5s'   '%5s'   '%8s'   '%10s'    '%10s'   '%10s' } ;
peas_varform        = {'%5.1f' '%4.1f' '%5.1f' '%0.5f' '%.1f' '%3.0f' '%.4d' '%3.0f'} ;
peas_bounds.xy(:,1) = [   10     .1    0.0001     0       0         0  0  0] ;
peas_bounds.xy(:,2) = [ 1000     .3    0.1     5000    1000     0.002 0  0] ;

% visco1d not yet implemented. Need to put proper names, etc
visco1d_names         = {'Len'  'Wid'    'Dep'   'Dip'   'Strike' 'xE'    'xN'    'ss'    'ds'    'op'   };   
visco1d_namform       = {'%5s'  '%5s'    '%5s'   '%6s'   '%7s'    '%5s'   '%5s'   '%5s'   '%5s'   '%5s'  };   
visco1d_varform       = {'%5.2f' '%5.2f' '%5.2f' '%6.2f' '%7.2f'  '%5.2f' '%5.2f' '%5.2f' '%5.2f' '%5.2f'};
visco1d_bounds.xy(:,1)= [  2       2       2     -89.9     0.1    x_start y_start  0       0       999   ];
visco1d_bounds.xy(:,2)= [ 10      10      10      89.9   179.9    x_end   y_end    0       0       999   ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create the concatenated strings and bounds vector %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[namform,varform,names,bounds]=deal([]);

%
% Dislocations
%

for i=1:double(N_disloc) 
   namform=[namform disloc_namform]; 
   varform=[varform disloc_varform]; 
   names  =[names   disloc_names  ];
end
for i=1:double(N_fault) 
   namform=[namform fault_namform]; 
   varform=[varform fault_varform]; 
   names  =[names   fault_names  ];
end
%
% Mogi,Penny,McTigue,Yang
%

for i=1:double(N_mogi)
   namform=[namform mogi_namform]; 
   varform=[varform mogi_varform]; 
   names  =[names   mogi_names  ];
end

if N_penny >=1
   namform=[namform penny_namform]; 
   varform=[varform penny_varform]; 
   names  =[names   penny_names  ];
end

if N_mctigue >=1
   namform=[namform mctigue_namform]; 
   varform=[varform mctigue_varform]; 
   names  =[names   mctigue_names  ];
end

if N_yang >= 1
   namform=[namform yang_namform]; 
   varform=[varform yang_varform]; 
   names  =[names   yang_names  ];
end

if N_multidisloc >= 1
   namform=[namform multidisloc_namform]; 
   varform=[varform multidisloc_varform]; 
   names  =[names   multidisloc_names  ];
end

if N_squaredisloc >= 1
   namform=[namform squaredisloc_namform]; 
   varform=[varform squaredisloc_varform]; 
   names  =[names   squaredisloc_names  ];
end

if N_visco1d >= 1
   namform=[namform visco1d_namform]; 
   varform=[varform visco1d_varform]; 
   names  =[names   visco1d_names  ];
end

if N_lockedandcreep >= 1
   namform=[namform lockedandcreep_namform]; 
   varform=[varform lockedandcreep_varform]; 
   names  =[names   lockedandcreep_names  ];
end

if N_peas >= 1
   namform=[namform peas_namform]; 
   varform=[varform peas_varform]; 
   names  =[names   peas_names  ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelopt.ParNames     = names;
modelopt.ParForm      = varform;
modelopt.ParNamesForm = namform;
modelopt.N_sources    = N_sources;

if N_disloc          default_bounds.disloc          = disloc_bounds;         end
if N_fault           default_bounds.fault           = fault_bounds;          end
if N_mogi            default_bounds.mogi            = mogi_bounds;           end
if N_penny           default_bounds.penny           = penny_bounds;          end
if N_mctigue         default_bounds.mctigue         = mctigue_bounds;        end
if N_yang            default_bounds.yang            = yang_bounds;           end
if N_multidisloc     default_bounds.multidisloc     = multidisloc_bounds;    end
if N_visco1d         default_bounds.visco1d         = visco1d_bounds;        end
if N_lockedandcreep  default_bounds.lockedandcreep  = lockedandcreep_bounds; end
if N_peas            default_bounds.peas            = peas_bounds          ; end
if N_squaredisloc    default_bounds.squaredisloc    = squaredisloc_bounds;   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N_visco1d  modelopt.visco1dopt.basemap  = basemap; end     %this is a temporay hack so that visco1d can generate lola coordinates
% xy to lola conversion for model par if given
if exist('par','var')
   if ~isfield(par,'lola') && isfield(par,'xy')   modelopt.par.lola=modelpar_lola2xy(par.xy,basemap,modelopt,-1);  end
   if ~isfield(par,'xy') && isfield(par,'lola')   modelopt.par.xy=modelpar_lola2xy(par.lola,basemap,modelopt, 1);  end
end

%TODO: here we clean out unnecessary fields
%if ~N_disloc modelopt=rmfield(modelopt,'dislocopt');  end      %%dislocopt.nu is used in ForwardModel
