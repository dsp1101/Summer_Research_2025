function [stop, options, optchanged] = trackProgress(options, optimvalues, flag, history)
% Output function for simulated annealing to track optimization progress
%
% This function is called by simulannealbnd at each iteration to record
% the misfit and parameter values over time
%
% Inputs:
%   options - Current options structure
%   optimvalues - Structure containing current optimization values
%   flag - Current state of the algorithm
%   history - OptimizationHistory object (handle class)
%
% Outputs:
%   stop - Boolean to stop optimization (always false here)
%   options - Modified options (not changed here)
%   optchanged - Boolean indicating if options changed (always false)

    stop = false;
    optchanged = false;
    
    if strcmp(flag, 'iter')
        % Add current iteration data to history
        history.addIteration(optimvalues.fval, optimvalues.x);
        
        
        if mod(history.iteration, 100) == 0
            fprintf('Iteration %d: Misfit = %.6f Time: %.4f minutes\n', ...
                history.iteration, optimvalues.fval, toc/60);
        end
    end
end