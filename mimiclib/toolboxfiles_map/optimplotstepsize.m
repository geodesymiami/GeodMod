function stop = optimplotstepsize(x,optimValues,state)
% OPTIMPLOTSTEPSIZE Plot step size at each iteration.
%
%   STOP = OPTIMPLOTSTEPSIZE(X,OPTIMVALUES,STATE) plots 
%   OPTIMVALUES.stepsize.
%
%   Example:
%   Create an options structure that will use OPTIMPLOTSTEPSIZE as the plot
%   function
%     options = optimset('PlotFcns',@optimplotstepsize);
%
%   Pass the options into an optimization problem to view the plot
%      fmincon(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[0 0],[],[],options)

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:46:33 $

persistent plotavailable
stop = false;

switch state
    case 'init'
        if isfield(optimValues,'stepsize')
            plotavailable = true;
        else
            plotavailable = false;
            title('Step Size: not available','interp','none');
        end
    case 'iter'
        if plotavailable
            if optimValues.iteration == 1
                % The 'iter' case is  called during the zeroth iteration, but
                % stepsize is still empty.  Start plotting at the first
                % iteration.
                plotstepsize = plot(optimValues.iteration,optimValues.stepsize,'kd', ...
                    'MarkerFaceColor',[1 0 1]);
                title(['Step Size: ',num2str(optimValues.stepsize)],'interp','none');
                xlabel('Iteration','interp','none');
                ylabel('Step size','interp','none');
                set(plotstepsize,'Tag','optimplotstepsize');
            else
                plotstepsize = findobj(get(gca,'Children'),'Tag','optimplotstepsize');
                newX = [get(plotstepsize,'Xdata') optimValues.iteration];
                newY = [get(plotstepsize,'Ydata') optimValues.stepsize];
                set(plotstepsize,'Xdata',newX, 'Ydata',newY);
                set(get(gca,'Title'),'String',['Step Size: ',num2str(optimValues.stepsize)]);
            end
        end
end
