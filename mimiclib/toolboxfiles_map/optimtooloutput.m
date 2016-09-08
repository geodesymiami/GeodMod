function stop = optimtooloutput(xOutputfcn,optimval,state,varargin);
%OPTIMTOOLOUTPUT OutputFcn used to interact between the GUI and solvers.
%  It is used to capture and react to changes in the state of the GUI 
%  (STOP, PAUSE, RESUME, KILL). It also manages warnings by displaying them
%  in Status and Results panel. It also sets the iteration number.
%
%   Private to OPTIMTOOL.

%   Copyright 2005-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/12/10 21:50:03 $

%Initialize
stop = false;
drawnow;
% Three buttons in the GUI
STOP = 0;
RUN_RESUME = 1;
PAUSE = 2;

optimtoolGui = com.mathworks.toolbox.optim.OptimGUI.getOptimGUI; % Get a handle to the GUI
if isempty(optimtoolGui)
    stop = true;
    return;
end

switch state
    case 'init'
        [msg,id] = lastwarn;
        if ~isempty(msg)
            % Make sure that exactly one newline character is at the end of 'msg'
            newLineIndex = regexp(msg,'\n');
            if isempty(newLineIndex) || newLineIndex(end) ~=length(msg)
                msg = sprintf('%s\n',msg);
            end
            optimtoolGui.appendResults(['Warning: ',msg]); % Append to the 'Status and Results' panel
            msg(regexp(msg,'\n')) = ' '; % Remove all newline characters before showing 'msg' in the dialog
            warndlg(msg,'Optimization Tool'); % Show a warning dialog box
        end
        % To avoid displaying warnings with same ID multiple times we store the warning IDs in appdata.
        % We show it only once in the 'Status and Results' panel.
        setappdata(0,'last_warning_id_for_optimtool',id);
        return; % Nothing to do now
    case {'iter', 'interrupt'}
        [msg,id] = lastwarn;
        % In addition to non-empty message the warning ID must be new (not in the appdata)
        if ~isempty(msg) && ~strcmp(id,getappdata(0,'last_warning_id_for_optimtool'))
        % Add a newline to the message if it does not have already
        newLineIndex = regexp(msg,'\n');
        if isempty(newLineIndex) || newLineIndex(end) ~=length(msg)
            msg = sprintf('%s\n',msg);
        end
        optimtoolGui.appendResults(['Warning: ',msg]);
        lastwarn('');
        end
        optimtoolGui.setIteration(value2RHS(optimval.iteration)); % Update iteration number in the GUI
        % Action based on run mode of GUI
        RunMode = optimtoolGui.getRunMode;
        switch RunMode
            case RUN_RESUME

            case STOP
                stop = true;
                return;
            case PAUSE
                fprintf('%s\n%s\n','OPTIMTOOL is paused. MATLAB Command prompt will', ...
                    'not be accessible until the optimization solver is completed.');
                % If in pause state keeping looping here.
                while true
                    drawnow
                    if isempty(com.mathworks.toolbox.optim.OptimGUI.getOptimGUI)
                        stop = true; % The GUI was closed by user
                        return;
                    end
                    mode = optimtoolGui.getRunMode;
                    if mode == STOP
                        stop = true;
                        return;
                    elseif mode == RUN_RESUME
                        break;
                    end
                end % End while
            otherwise % Safegaurd
                return;
        end
    case 'done'
        optimtoolGui.setIteration(value2RHS(optimval.iteration)); % Update iteration number in the GUI
        if isappdata(0,'last_warning_id_for_optimtool')
            rmappdata(0,'last_warning_id_for_optimtool');
        end
        return;
    otherwise % Safegaurd
        return;
end
