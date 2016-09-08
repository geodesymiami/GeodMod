function [err,x,fval,exitMessage,nrow,ncolx,ncolf] = optimguirun(hashProb, hashOpt)
%OPTIMGUIRUN Optimization Toolbox GUI 'Start' button callback function.
%   Two arguments 'hashProb' and 'hashOpt' are Java hash tables for 
%   problem model and options model respectively.
%   The output 'err' is the error string returned by either readOptimHashTable
%   or callSolver functions. x,fval,exitMessage are outputs from the solver
%   and [nrow,ncol] is the size of 'x' which is needed to display 'x' in the GUI. 

%   Copyright 2005-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2007/12/10 21:50:40 $

err = '';
x = '';
fval = '';
exitMessage = '';
% Size of the result vector 'X' and 'fval' (nrows is always the number of
% solution)
nrow = [];
ncolx = [];
ncolf = [];
% Maximum length of 'X' vector to be shown in the GUI
MAX_NUM_ELEMENT_SHOW = 100;
% Get modified fields from the GUI
[probStruct,optStruct,errProb,errOpt] = readOptimHashTable(hashProb, hashOpt);
if ~isempty(errProb)
    err = errProb;
    return;
elseif ~isempty(errOpt)
    err = errOpt;
    return;
end
% Add appropriate output function for the 'solver' so
% that the GUI and the solver can interact in each iteration
optStruct = addOptimguiOutputFcn(optStruct,probStruct.solver);
% Store current warning state
[lastmsg, lastid] = lastwarn;
lastwarn('');

try
    % Call solver and save the result structure to MATLAB workspace (appdata)
    resultStruct = callSolver(probStruct,optStruct);
    setappdata(0,'optimTool_results_121677',resultStruct);
    exitMessage = resultStruct.output.message;
    x = resultStruct.x;
    
    % Set iteration number in the GUI
    optimGUI = com.mathworks.toolbox.optim.OptimGUI.getOptimGUI; % Get a handle to the GUI
    if ~isempty(optimGUI)
        % Update iteration number in the GUI
        try
            optimGUI.setIteration(value2RHS(resultStruct.output.iterations));
        catch % GA solvers have 'generations' and not 'iterations'
            optimGUI.setIteration(value2RHS(resultStruct.output.generations));
        end
    end

    % Least square solvers have 'resnorm' and not 'fval'
    try
        fval = resultStruct.fval;
    catch
        fval = resultStruct.resnorm;
    end
    % Return argument fval
    if ndims(fval) < 3 && numel(fval) <= MAX_NUM_ELEMENT_SHOW
        if isscalar(fval) % Single objective
            ncolf = 0; % GUI expects one less
        else % Multiple solution (only gamultiobj)
          [unused1, ncolf] = size(fval);
        end
    else
        ncolf = -1;
        fval  = [];
    end
    % Return argument x
    if ndims(x) < 3 && (isnumeric(x) || isa(x,'double')) && ...
            numel(x) <= MAX_NUM_ELEMENT_SHOW
            [nrow, ncolx] = size(x);          
    else
        nrow  = -1;
        ncolx = -1;
        x = [];
    end
    % Save the random states if it is returned in the resultStruct
    if isfield(resultStruct.output,'randstate')
        probStruct.randstate = resultStruct.output.randstate;
        probStruct.randnstate = resultStruct.output.randnstate;
        setoptimrandstates(probStruct,java.util.Hashtable,false);
    end
catch
    nrow = -1;
    ncolx = -1;
    ncolf = -1;
    x = [];
    err = lasterr;
end

% Restore the warning if there were no warnings from the solver
if isempty(lastwarn)
    lastwarn(lastmsg,lastid);
end

