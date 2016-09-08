function err = optimguiGeneratemfile(hashProb,hashOpt)
%optimguiGeneratemfile generates an M-file from OPTIMTOOL.
%   hashProb and hashOpt are Java hash tables containing information about
%   the problem and options model. hashProb and hashOpt contain only information 
%   that user has changed since last time the data (Java model) from the GUI 
%   was passed to MATLAB workspace. (E.g. at the time of exporting, running,
%   generating code.) This function will update the MATLAB workspace and will 
%   call generateMfile.m to generate an M-file.
%
%   Private to OPTIMTOOL

%   Copyright 2005-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/11/09 20:57:22 $

err = '';
% Get modified fields from the GUI
[probStruct,optStruct,errProb,errOpt] = readOptimHashTable(hashProb, hashOpt);
if ~isempty(errProb)
    err = errProb;
    return;
elseif ~isempty(errOpt)
    err = errOpt;
    return;
end

% Generate M-file with for modified problem and options structure
generateMfile(probStruct,optStruct);

%----------------------------------------------------------------------
function msg = generateMfile(probStruct,optStruct)
%generateMfile helper function for optimguiGeneratemfile. 
%
%   This function takes a problem structure 'probStruct' and an options
%   structure 'optStruct' and generates an M-file. The generated code sets
%   the non-default options and also calls the appropriate solver.
%

msg  = '';
% headerCode is the code for the generated function signature
headerCode = 'function [';
% optionsCode is the code for options settings
optionsCode = '';
% solverCode is the code for calling the appropriate solver
solverCode = '[';

% Get file name to use, remember the directory name
filespec = '*.m';
[optimMfileName,pn] = uiputfile(filespec,'Generate M-File','untitled.m');
if isequal(optimMfileName,0) || isequal(pn,0)
    return
end
if ~ismember('.',optimMfileName)
    optimMfileName = [optimMfileName '.m'];
end
optimMfileName = sprintf('%s%s',pn,optimMfileName);


% Get M file name with .m suffix, and get corresponding function name
if length(optimMfileName)<2 || ~isequal(optimMfileName(end-1:end),'.m')
    optimMfileName = sprintf('%s.m',optimMfileName);
end
[dirname,fcnname] = fileparts(optimMfileName);

solver = probStruct.solver; 
% Problem structure will tell us about the number of input arguments and
% also about the order in which the solver accept them. The input 'probStruct'
% contains fields for all the solvers. The output 'probStruct' will only 
% contain fields relevant to 'solver'.
probStruct = createProblemStruct(solver,[],probStruct);
% We do not want the 'solver' field in the problem structure
probStruct = rmfield(probStruct,'solver');
% Check if random states are used in M-file generation
if isfield(probStruct,'randstate')
    probStruct = rmfield(probStruct,{'randstate','randnstate'});
end

problemFields = fieldnames(probStruct);
problemNumFields = length(problemFields);

% Result structure will tell us about the number of output arguments and
% also about the order in which it comes from the solver.
resultStruct = createResultsStruct(solver);
resultsFields = fieldnames(resultStruct);
resultsNumFields = length(resultsFields);
% Append solver's output list to the  signature of the generated 
% function and also to the calling syntax for the solver.
for i = 1:resultsNumFields
   headerCode = [headerCode sprintf('%s,',resultsFields{i})]; 
   solverCode = [solverCode sprintf('%s,',resultsFields{i})]; 
end
% Replace the last char ',' (comma) by ']' in the headerCode/solverCode 
headerCode(end) = ']'; solverCode(end) = ']';
% Add '=' character to the end
headerCode = [headerCode ' = ']; solverCode = [solverCode ' = '];
% Add the file name to headerCode
headerCode = [headerCode sprintf('%s(',fcnname)];
% Add the 'solver' name to solverCode
solverCode = [solverCode sprintf('%s\n','...') sprintf('%s(',solver)];

% For each field of the problem structure
for i = 1:problemNumFields
    probField = problemFields{i};
    probValue = probStruct.(probField);
    tempcode = getStringValue(probValue,probField,true);
    solverCode = [solverCode tempcode sprintf(',')];
        % If value is numeric and non-empty then this code also goes into
        % the headerCode.
        if ~isempty(probValue)  && isnumeric(probValue)
            headerCode = [headerCode tempcode sprintf(',')];  
        end
end

% After all the problem data is written we want to add options as the last
% argument to solverCode
solverCode = [solverCode sprintf('options);')];
% Get fieldnames of the options structrue
[optStruct,optionsFcn] = createOptionsStruct(solver,optStruct);
optionsFields = fieldnames(optStruct);
optionsNumFields = length(optionsFields);

% Start with default options
optionsCode = [optionsCode sprintf('\n%s\n','% Start with the default options')];
optionsCode = [optionsCode sprintf('options = %s;\n',optionsFcn)];

% For each property
optionsCode = [optionsCode sprintf('%s\n','% Modify options setting')];
for i = 1:optionsNumFields
    optField = optionsFields{i};
    optValue = optStruct.(optField);
    if  ~isempty(optValue)  % don't generate code for defaults.
        tempcode = getStringValue(optValue,optField,false);
        optionsCode = [optionsCode sprintf('options = %s(options,''%s'' ,%s);\n', ...
            optionsFcn,optField,tempcode)];
        % If value is numeric and non-empty then this code also goes into
        % the headerCode.
        if ~isempty(optValue) && isnumeric(optValue)
            headerCode = [headerCode tempcode sprintf(',')];    
        end
    end
end

% Replace the last char ',' (comma) by the ')' (closing bracket) in the headerCode
if headerCode(end) == ','
    headerCode(end) = ')';
else
    headerCode(end) = '';
end
% Add some generic comments in the file
headerCode = [headerCode sprintf('\n%s\n','% This is an auto generated M-file from Optimization Tool.')];
%headerCode = [headerCode sprintf('%s\n','% Optimization Toolbox is required to run this M-file.')];

% Concatenate all the generated code into one
code = [sprintf('%s',headerCode),sprintf('%s',optionsCode),sprintf('%s',solverCode)];
% Open and write to the file
[fid,message] = fopen(optimMfileName,'w');
if fid==-1
    msg = sprintf('Error trying to write to %s:\n%s',optimMfileName,message);
    errordlg(msg,'Error Saving M File','modal');
    return
end
% Close the file
fprintf(fid,'%s\n',code);
st = fclose(fid);
if st ~= 0
    msg = sprintf('%s%s','Error closing file ',fcnname);
    return;
end
% Open the M file just created
edit(optimMfileName)

function rhsValue = getStringValue(value,FieldName,problem_type)
% Function to wrap around value2RHS so that matrices are not
% converted to strings. 
% If problem_type is 'true' then string for solver's input  arguments are generated
% otherwise for options using the 'FieldName'

if ~isempty(value) && isnumeric(value)
    if problem_type
        rhsValue = sprintf('%s',FieldName); % This is for problem values that are non-empty
    else
        rhsValue = sprintf('%s_Data',FieldName); % This is for options values that are matrices/vectors
    end
else
    rhsValue = value2RHS(value);
end

