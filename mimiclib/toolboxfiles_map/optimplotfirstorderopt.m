function stop = optimplotfirstorderopt(x,optimValues,state)
% OPTIMPLOTFIRSTORDEROPT Plot first-order optimality at each iteration.
%
%   STOP = OPTIMPLOTFIRSTORDEROPT(X,OPTIMVALUES,STATE) plots
%   OPTIMVALUES.firstorderopt.
%
%   Example:
%   Create an options structure that will use OPTIMPLOTFIRSTORDEROPT as the
%   plot function
%     options = optimset('PlotFcns',@optimplotfirstorderopt);
%
%   Pass the options into an optimization problem to view the plot
%      fmincon(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[0 0],[],[],options)

%   Copyright 2006-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/12/10 21:50:01 $

persistent plotavailable
stop = false;

switch state
    case 'iter'
        if optimValues.iteration == 1
            if isfield(optimValues,'firstorderopt') && ~isempty(optimValues.firstorderopt)
                plotavailable = true;

                % The 'iter' case is  called during the zeroth iteration, but
                % firstorderopt may still  be empty.  Start plotting at the
                % first iteration.
                plotfirstorderopt = plot(optimValues.iteration,optimValues.firstorderopt,'kd', ...
                    'MarkerFaceColor',[1 0 1]);
                title(['First-order Optimality: ',num2str(optimValues.firstorderopt)],'interp','none');
                xlabel('Iteration','interp','none');
                ylabel('First-order optimality','interp','none');
                set(plotfirstorderopt,'Tag','optimplotfirstorderopt');
            else % firstorderopt field does not exist or is empty
                plotavailable = false;
                title('First-order Optimality: not available','interp','none');
            end
        else
            if plotavailable
                plotfirstorderopt = findobj(get(gca,'Children'),'Tag','optimplotfirstorderopt');
                newX = [get(plotfirstorderopt,'Xdata') optimValues.iteration];
                newY = [get(plotfirstorderopt,'Ydata') optimValues.firstorderopt];
                set(plotfirstorderopt,'Xdata',newX, 'Ydata',newY);
                set(get(gca,'Title'),'String',['First-order Optimality: ',num2str(optimValues.firstorderopt)]);
            end
        end
end