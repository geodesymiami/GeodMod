function [selection, problemModel, optionModel] = optimguiImportProblem()
%optimguiImportProblem Optimtool helper function to import problem structure. 
% It presents a list dialog to the user. The dialog contains list of valid optim 
% problem structures. This function saves the problem and options structure (which 
% is a field of problem structure) to the MATLAB workspace and also returns the 
% equivalent Java hash tables 'problemModel' and 'optionModel' to the GUI. The name 
% of the selected variable 'selection' to be imported is also returned to the GUI.
%
%   Private to OPTIMTOOL

%   Copyright 2005-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $  $Date: 2007/12/10 21:50:38 $


selection = '';
optionModel = '';
problemModel = '';
names = {};
optionsFieldnames = fieldnames(createOptionsStruct('all'));
% We check for one matching field (to make sure that all optimization
% solvers get through this test)
minNumberOfOptionsToCheck = 1; % Display

probFieldnames = fieldnames(createProblemStruct('all',[]));
% Required field names for problem structure
requiredFields = {'solver','options'};
validValues    =  {fieldnames(createProblemStruct('solvers')), {} };

whoslist =evalin('base','whos');
for i = 1:length(whoslist)
    if strcmp(whoslist(i).class, 'struct') && strcmp(num2str(whoslist(i).size), '1  1')
        s = evalin('base', whoslist(i).name);
        if validOptimProblemStruct(s,requiredFields,validValues) && ...
                validOptions(s.options,optionsFieldnames,minNumberOfOptionsToCheck)
            names{end + 1} = whoslist(i).name;
        end
    end
end
   
if isempty(names) 
    msgbox('There are no problem structures in the workspace.', 'Optimization Tool');
else
    [selectionIndex, Answer] = listdlg('ListString', names, 'SelectionMode', 'Single', ...
            'ListSize', [250 200], 'Name', 'Import Optimization Problem', ...
            'PromptString', 'Select a problem structure to import:', ...
            'OKString', 'Import');
    % Answer == 1 means that user pressed the 'Import' button
    if Answer == 1
        selection = names{selectionIndex};
        [unused1,unused2,probStruct] = validOptimProblemStruct(evalin('base', selection),requiredFields,validValues);
        options = createOptionsStruct('all',probStruct.options); 
        probStruct = rmfield(probStruct,'options');
        probStruct = createProblemStruct('all',[],probStruct);
        % Stuff all the fields into the hashtable.
        problemModel = createHashTable(probStruct,probFieldnames,selection); % Create Java hashtable
        optionModel = createHashTable(options,optionsFieldnames,[selection,'.options']);
        % Save problem and options to the MATLAB workspace (appdata)
        setappdata(0,'optimTool_Problem_Data',probStruct);
        setappdata(0,'optimTool_Options_Data',options);
     end
end    

% Reset Java hashtable for options and problem change
resetOptimtoolHashTable('optimTool_Problem_HashTable');
resetOptimtoolHashTable('optimTool_Options_HashTable');
