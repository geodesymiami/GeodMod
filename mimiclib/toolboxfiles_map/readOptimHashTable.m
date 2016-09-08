function [probStruct,optStruct,errProb,errOpt] = readOptimHashTable(hashProb, hashOptions)
%readOptimHashTable Read hash table from optimtool and return MATLAB structures.
%  Private to OPTIMTOOL

%   Copyright 2005-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/12/10 21:50:06 $

% Add problem/options keys that are workspace variables and have been used
% previously. Call to this helper function will ensure that changes in
% workspace variables are always picked up.
[hashProb,hashOptions] = addModifiableKeys(hashProb,hashOptions);

% Read only the modified options keys (so far!) from the hash table
modifiedOptionsKeys = toArray(keySet(hashOptions));

% Get the saved problem structure from MATLAB workspace (appdata).
probStruct = getappdata(0,'optimTool_Problem_Data');
errProb = '';
% We get values of all the problem keys for the selected 'solver' from the
% hash table. The 'solver' is present first time and every time the
% 'solver' field changes.
currentSolver = hashProb.get('solver'); 

% Get valid fields for the currentSolver
currentprobFields = fieldnames(createProblemStruct(currentSolver));
for i = 1:length(currentprobFields)
    [probStruct, errs] = getData(hashProb,currentprobFields{i},probStruct);
    % Concatenate errors from all keys into one to be shown in the GUI just once.
    errProb = [errProb errs];
end

% Special code to handle states for random number generators (problem
% structures only). problemStruct gets the random states for the
% appropriate solver from the saved appdata
if any(strcmpi(probStruct.solver,{'ga','patternsearch','simulannealbnd', ...
        'threshacceptbnd','gamultiobj'}))
    [probStruct,errs] = getoptimrandstates(probStruct,hashProb);
    errProb = [errProb errs];
end

optStruct = getappdata(0,'optimTool_Options_Data');
errOpt = '';
% We get values for modified keys of the options hash table in a loop.
for i = 1:length(modifiedOptionsKeys)
    [optStruct, errs] = getData(hashOptions,modifiedOptionsKeys(i),optStruct);
    % We do not care about errors coming from disabled options
    if ~isempty(errs)
        optimGUI = com.mathworks.toolbox.optim.OptimGUI.getOptimGUI; % Get a handle to the GUI
        if ~isempty(optimGUI) && ~optimGUI.isOptionEnabled(modifiedOptionsKeys(i))
            % Get the old value which is stored in the appdata
            tempStruct = getappdata(0,'optimTool_Options_Data');
            optStruct.(modifiedOptionsKeys(i)) = tempStruct.(modifiedOptionsKeys(i));
            continue;
        end
    end
    % Concatenate errors from all keys into one to be shown in the GUI just once.
    errOpt = [errOpt errs];  
end
% Update problem and options stored in the MATLAB workspace (appdata)
setappdata(0,'optimTool_Problem_Data',probStruct);
setappdata(0,'optimTool_Options_Data',optStruct);
%---------------------------------------------------------
% Function to add all modifiable keys to hash table before their values are 
% evaluated in workspace
function [problemHash,optionsHash] = addModifiableKeys(hashProb,hashOptions)

% Check if problem hash table is already stored
if isappdata(0,'optimTool_Problem_HashTable')
    problemHash = getappdata(0,'optimTool_Problem_HashTable');
else
    problemHash = hashProb;
end

% Check if options hash table is already stored
if isappdata(0,'optimTool_Options_HashTable')
    optionsHash = getappdata(0,'optimTool_Options_HashTable');
else
    optionsHash = hashOptions;
end

newProbKeys = toArray(keySet(hashProb));
newOptsKeys = toArray(keySet(hashOptions));
% We get values for modified keys of the problem hash table in a loop.
for i = 1:length(newProbKeys)
  problemHash.put(newProbKeys(i),hashProb.get(newProbKeys(i)));
end
% We get values for modified keys of the options hash table in a loop.
for i = 1:length(newOptsKeys)
  optionsHash.put(newOptsKeys(i),hashOptions.get(newOptsKeys(i)));
end
% Update problem and options stored in the MATLAB workspace (appdata)
setappdata(0,'optimTool_Problem_HashTable',problemHash);
setappdata(0,'optimTool_Options_HashTable',optionsHash);

%---------------------------------------------------------
% getData function to interpret the values in hashtable
function [myStruct, err] = getData(myHash,key,myStruct)

err = '';
if myHash.containsKey(key)
    value = myHash.get(key);
else
    return;
end

type = dataTypeOfKey(key); % What are the types of value the key can take

% Need to evaluate the value if not knwon (other than []).
% If 'type' is 'string' then it can not be evaluated in MATLAB base (hence try-catch).
% If 'type' is 'stringAndValue' then need to evaluate the when it is a value e.g., cell array.
if strcmp(value,'[]') || ~strcmp(type,'string')
    try
        value = evalin('base', value);
    catch 
        % if the type is 'stringAndValue' then evalin will error when it is
        % a string so keep the string as it is
        if strcmp(type,'stringAndValue')
            myStruct.(key) = value;
            return;
        end
        % Check if the value is already in the structure (encrypted form)
        if ~isempty([strfind(value,'<userStructure>') ...
                strfind(value,'<userClass>') ...
                strfind(value,'<userData>')]);
            return;
        end
        err = sprintf('Error in %s: %s \n', key, lasterr);
    end
end
% Update the structure with the field and the value
myStruct.(key) = value;

%---------------------------------------------
% isValueString check if the optName can be a string
function type = dataTypeOfKey(optName)
% This function returns the type of field 'optName'. The possible types are 'string', 
% 'mixed', and 'other'.

% The possible values of these fields are strings only.
stringOnly =   {'Simplex','NonlEqnAlgorithm', ...
    'NodeSearchStrategy','MeritFunction', ...
    'LineSearchType','LevenbergMarquardt', ...
    'LargeScale','Jacobian', ...
    'InitialHessType','HessUpdate', ...
    'GradObj','GradConstr','Diagnostics', ...
    'DerivativeCheck','BranchStrategy', ...
    'FunValCheck','Display','Vectorized', ...
    'PopulationType','MigrationDirection', ...
    'MeshAccelerator','MeshRotate', ...
    'ScaleMesh','PollMethod','CompletePoll','PollingOrder', ...
    'CompleteSearch','Cache','DataType','solver', ...
    'AlwaysHonorConstraints','FinDiffType','Algorithm', ...
    'ScaleProblem','SubproblemAlgorithm'};

% This group of options (mixed) is currently not used in the function 'getData'.
% The possible values of these fields may contain either fixed strings or numbers.
stringAndValue = {'TypicalX','MaxSQPIter', ...
    'MaxRLPIter','MaxPCGIter', ...
    'MaxNodes','JacobPattern', ...
    'HessPattern','MaxIter', ...
    'MaxFunEvals','Hessian','TolX'};

if ismember(optName,stringOnly)
    type = 'string';
elseif ismember(optName,stringAndValue)
    type = 'stringAndValue';
else 
    type = 'other';
end
