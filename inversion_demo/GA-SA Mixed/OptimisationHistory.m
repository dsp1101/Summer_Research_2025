classdef OptimisationHistory < handle
    % OptimizationHistory - Handle class to track optimization progress
    %
    % Stores history of misfit values and parameters during simulated
    % annealing optimisation
    % 
    % % Usage:
    %   hist = OptimizationHistory();
    %   hist.addIteration(misfit, params);
    %   misfits = hist.getMisfits();
    properties
        misfit
        params
        iteration = 0
        maxSize
        nParams = []
        growthFactor = 1.5
    end

    methods
        function obj=OptimisationHistory(estimatedIterations)
        % Input = estimatedIterations (optional) - rough guess of total
        % number of iterations
        if nargin<1
            estimatedIterations = 10000;
        end
        
        obj.maxSize = estimatedIterations;
        obj.misfit = zeros(estimatedIterations, 1);
        
        end

        function addIteration(obj, misfitValue, paramValues)
            % Add a new iteration's data to the history
            if obj.iteration == 0
                obj.nParams=length(paramValues);
                obj.params=zeros(obj.maxSize, obj.nParams);
            end

            obj.iteration = obj.iteration + 1;

            if obj.iteration>obj.maxSize
                obj.growArrays();
            end
            
            obj.misfit(obj.iteration) = misfitValue;
            obj.params(obj.iteration, :) = paramValues;
        end

        function growArrays(obj)
            newSize=ceil(obj.maxSize*obj.growthFactor);

            obj.misfit(newSize)=0;
            obj.params(newSize, :) = 0;
            obj.maxSize = newSize;
        end

        function m = getMisfits(obj)
            m = obj.misfit(1:obj.iteration);
        end

        function p = getParams(obj)
            % Get all parameter values
            p = obj.params(1:obj.iteration,:);
        end
        
        function n = getIterationCount(obj)
            % Get total number of iterations
            n = obj.iteration;
        end
        
        function reset(obj, estimatedIterations)
            % Reset the history
            if nargin<2
                estimatedIterations = 10000;
            end
            obj.iteration=0;
            obj.maxSize
            obj.misfit = estimatedIterations;
            obj.misfit = zeros(estimatedIterations, 1);
            obj.params = zeros(estimatedIterations, 6);
        end

    end
end