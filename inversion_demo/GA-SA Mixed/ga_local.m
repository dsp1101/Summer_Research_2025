%% Helper Function to perform GA optimisation within SA GA mixed function%%
%
% GA_local-Genetic Algorithm implementation which keeps initial population
%          within a spread of given initial paramaters
%
%
%
%
% Written by Daniel Spencer
% 23/01/26
% dsp65@uclive.ac.nz // daniel.spencer2007@gmail.com

function [population, misfit, simSpec, f, bestIndiv, globalBestIndiv, ...
        finalGen] = ga_local(popSize, nGen, mutationRate, crossoverRate, ...
        eliteCount, tournamentSize, geomParams, geomFlag, srcParams, ...
        srcFlag, srcStyle, lowerBnds, upperBnds, discrParams, temp, ...
            filterProps, data_freq, freqLim, useTimeLimit, maxTime, initialParams)

% INPUTS:
%   popSize         - population size (number of individuals)
%   nGen            - number of generations
%   mutationRate    - probability of mutation (0-1)
%   crossoverRate   - probability of crossover (0-1)
%   eliteCount      - number of elite individuals to preserve
%   tournamentSize  - size of tournament for selection
%   geomParams      - initial geometry parameters [depth, r1, r2, ...]
%   geomFlag        - boolean, invert for geometry (1) or not (0)
%   srcParams       - source parameters
%   srcFlag         - boolean, invert for source (1) or not (0)
%   srcStyle        - source model style (e.g., 'Brune')
%   lowerBnds       - lower bounds for parameters
%   upperBnds       - upper bounds for parameters
%   discrParams     - discretization parameters [T N Nf Nyquist dt]
%   temp            - temperature [crater, atmosphere]
%   filterProps     - filter properties [band, order, Fs]
%   data_freq       - observed data [frequency, amplitude]
%   freqLim         - high frequency limit for misfit calculation
%   useTimeLimit    - boolean, use time limit (1) or generation limit (0)
%   maxTime         - maximum time in seconds (if useTimeLimit=1)
%
% OUTPUTS:
%   population      - final population (normalized parameters in [0,1])
%   misfit          - best misfit per generation
%   simSpec         - simulated spectrum from best individual
%   f               - frequency vector
%   bestIndiv       - best individual per generation (normalized)
%   globalBestIndiv - best individual across all generations
%   finalGen        - number of completed generations

%% Determine number of parameters
nParamsGeom = geomFlag * length(geomParams);
nParamsSrc  = srcFlag * length(srcParams);
nParams     = nParamsGeom + nParamsSrc;

%% Initialize population

%  Validate initial parameters are within bounds
initialParams = max(lowerBnds, min(upperBnds, initialParams));

% normalize initial parameters to [0,1]
initialParamsNorm = (initialParams - lowerBnds) ./ (upperBnds - lowerBnds);

%  2: Increase spread slightly for better diversity
spread = 0.05;  % Increased from 0.02
population = repmat(initialParamsNorm, popSize, 1) + spread * randn(popSize, nParams);

% enforce [0,1] bounds
population = max(0, min(1, population));

% CRITICAL: Ensure first individual is EXACTLY the initial parameters
% This guarantees SA's best solution is preserved in the population
population(1, :) = initialParamsNorm;

%  3: Ensure initial population satisfies monotonic decreasing constraint
% Skip individual 1 since it's the validated initial parameters
if geomFlag && nParamsGeom > 1
    for i = 2:popSize  % Start from 2, skip the exact initial params at index 1
        % Get physical radii for this individual
        params_phys = lowerBnds + population(i,:) .* (upperBnds - lowerBnds);
        radii = params_phys(2:nParamsGeom);
        
        % If not monotonically decreasing, sort them
        if any(diff(radii) >= 0)
            radii_sorted = sort(radii, 'descend');
            % Add small decreasing perturbations to ensure strict decrease
            for j = 2:length(radii_sorted)
                if radii_sorted(j) >= radii_sorted(j-1)
                    radii_sorted(j) = radii_sorted(j-1) - 0.1;
                end
            end
            % Ensure within bounds
            radii_sorted = max(lowerBnds(2:nParamsGeom), min(upperBnds(2:nParamsGeom), radii_sorted));
            % Convert back to normalized
            population(i, 2:nParamsGeom) = (radii_sorted - lowerBnds(2:nParamsGeom)) ./ ...
                                            (upperBnds(2:nParamsGeom) - lowerBnds(2:nParamsGeom));
        end
    end
end

misfit = zeros(nGen,1);
bestIndiv = zeros(nGen, nParams);
fitness = zeros(popSize,1);

%% Extract discretization & temperature parameters
T = discrParams(1);
N = discrParams(2);
Nf = discrParams(3);
Nyquist = discrParams(4);
dt = discrParams(5);
freq = [0 Nyquist];

craterTemp = temp(1);
atmoTemp = temp(2);
M = problemParametersInv(craterTemp, atmoTemp);

%% Filter & data
filterband = filterProps(1:2);
filterorder = filterProps(3);
Fs = filterProps(4);

dataSpectrum.f = data_freq(:,1);
dataSpectrum.P = abs(data_freq(:,2))./max(abs(data_freq(:,2)));

% Pre-compute frequency index for band limiting
idx = find(dataSpectrum.f > freqLim, 1, 'first');
if isempty(idx)
    idx = length(dataSpectrum.f);
end

fprintf('Starting GA inversion: %d individuals, %d generations, %.0f max seconds...\n', popSize, nGen, maxTime);
fprintf('Initial parameters: ');
fprintf('%.2f ', initialParams);
fprintf('\n');

%% GA loop %%
tStart = tic;
gen = 0;
globalBestMisfit = Inf;  %  4: Initialize to Inf instead of []
globalBestIndiv = initialParamsNorm;

while true
    gen = gen + 1; 
    % --- Evaluate population ---
    for i = 1:popSize
        % Scale parameters to physical bounds
        params_scaled = lowerBnds + population(i,:) .* (upperBnds - lowerBnds);

        % Geometry and source extraction
        if geomFlag
            geomParams_current = params_scaled(1:nParamsGeom);
        else
            geomParams_current = geomParams;
        end
        if srcFlag
            srcParams_current = params_scaled(end);
        else
            srcParams_current = srcParams;
        end

        % Forward model evaluation
        try
            % Crater shape
            shape = geomFunction(geomParams_current);
            depth = shape(1,1);
            
            % Validate geometry
            if depth <= 0 || any(shape(:,2) <= 0)
                fitness(i) = 1e10; continue;
            end
            
            %  5: Relax constraint slightly - allow equal values with small tolerance
            if geomFlag && nParamsGeom > 1
                radii = geomParams_current(2:end);
                % Use <= instead of < to allow tiny differences
                if any(diff(radii) > 0.01)  % Relaxed from strictly < 0
                    fitness(i) = 1e10; continue;
                end
            end

            % Source function
            [S, ~, ~, ~] = sourceFunction(1, srcParams_current, srcStyle, discrParams);

            % Acoustic transfer
            style = 'baffled piston';
            order = 4;
            res = resonance1d(shape, depth, freq, Nf, style, order, M);

            % Convolve transfer function and source
            sim.f = res.f;
            sim.P = (res.P(1:N/2+1).*S(1:N/2+1))./max(res.P(1:N/2+1).*S(1:N/2+1));

            % Convert to time domain
            sim.pFreq = [sim.P conj(sim.P(end-1:-1:2))];
            sim.pFreq(N/2+1) = real(sim.pFreq(N/2+1));
            sim.pTime = ifft(sim.pFreq,'symmetric')/dt;

            % Apply filtering
            sim.pTimeFilt = bandpass_butterworth(sim.pTime, filterband, Fs, filterorder);

            % Convert back to frequency domain
            L = length(sim.pTimeFilt);
            NFFT = 2^nextpow2(L);
            sim.pFreqFilt = fft(sim.pTimeFilt, NFFT) / L;
            sim.pFreqAbsNorm = abs(sim.pFreqFilt(1:NFFT/2+1))./max(abs(sim.pFreqFilt(1:NFFT/2+1)));
            sim.fFilt = Fs/2*linspace(0,1,NFFT/2+1);

            % Interpolate simulation onto same frequency vector as data
            sim.pFreqIntp = pchip(sim.fFilt, sim.pFreqAbsNorm, dataSpectrum.f);
            
            % Validate interpolation
            if any(isnan(sim.pFreqIntp)) || any(isinf(sim.pFreqIntp))
                fitness(i) = 1e10; continue;
            end

            % Band limit signal
            dataAmpSpec = abs(dataSpectrum.P(1:idx));
            simAmpSpec = abs(sim.pFreqIntp(1:idx));

            % Use L2 norm
            fitness(i) = norm(dataAmpSpec - simAmpSpec);
            
            % Validate fitness
            if isnan(fitness(i)) || isinf(fitness(i))
                fitness(i) = 1e10;
            end

        catch ME
            fitness(i) = 1e10;
        end
    end

    %% Store best individual
    [currentBest, bestIdx] = min(fitness);
    
    if currentBest < globalBestMisfit  %  6: Simpler comparison
        globalBestMisfit = currentBest;
        globalBestIndiv = population(bestIdx,:);
    end
    
    misfit(gen) = currentBest;
    bestIndiv(gen,:) = population(bestIdx,:);

    %% Selection - tournament
    newPop = zeros(popSize, nParams);
    [~, sortedIdx] = sort(fitness);
    newPop(1:eliteCount,:) = population(sortedIdx(1:eliteCount),:);
    for j = eliteCount+1:popSize
        tIdx = randperm(popSize, tournamentSize);
        [~, wLocal] = min(fitness(tIdx));
        winner = tIdx(wLocal);
        newPop(j,:) = population(winner,:);
    end

    %% Crossover
    for j = eliteCount+1:2:popSize-1
        if rand < crossoverRate && nParams > 1
            cp = randi([1,nParams-1]);
            tmp = newPop(j,cp+1:end);
            newPop(j,cp+1:end) = newPop(j+1,cp+1:end);
            newPop(j+1,cp+1:end) = tmp;
        end
    end

    %% Mutation
    for j = eliteCount+1:popSize
        mutMask = rand(1,nParams) < mutationRate;
        newPop(j,mutMask) = rand(1,sum(mutMask));
    end
    
    population = newPop;
    
    %% Print progress
    if mod(gen,50)==0 || gen==1
        if useTimeLimit
            fprintf('Generation %d | Elapsed: %.2f min | Best Misfit: %.6f\n', ...
                    gen, toc(tStart)/60, misfit(gen));
        else
            fprintf('Generation %d/%d | Best Misfit: %.6f\n', gen, nGen, misfit(gen));
        end
        globalBestPhysical = lowerBnds + globalBestIndiv .* (upperBnds - lowerBnds);
        fprintf('  Best params: ');
        fprintf('%.2f ', globalBestPhysical);
        fprintf('\n');
    end
    
    %% Termination check
    if useTimeLimit
        if toc(tStart) >= maxTime
            finalGen = gen;
            break;
        end
    else
        if gen >= nGen
            finalGen = gen;
            break;
        end
    end
end

%% Trim output arrays to actual number of generations
misfit = misfit(1:finalGen);
bestIndiv = bestIndiv(1:finalGen,:);

%%  7: Convert final parameters to physical values
params_final = lowerBnds + globalBestIndiv .* (upperBnds - lowerBnds);

fprintf('\n--- GA Refinement Complete ---\n');
fprintf('Final best misfit: %.6f\n', globalBestMisfit);
fprintf('Generations completed: %d\n', finalGen);
fprintf("Initial Parameters: ")
fprintf('%.2f', initialParams)
fprintf('\n')
fprintf('Final parameters: ');
fprintf('%.2f', params_final);
fprintf('\n');

if geomFlag
    geomParams_final = params_final(1:nParamsGeom);
else
    geomParams_final = geomParams;
end
if srcFlag
    srcParams_final = params_final(end);
else
    srcParams_final = srcParams;
end

% Forward model for best individual
shape = geomFunction(geomParams_final);
depth = shape(1,1);
[S, ~, ~, ~] = sourceFunction(1, srcParams_final, srcStyle, discrParams);

style = 'baffled piston';
order = 4;
res = resonance1d(shape, depth, freq, Nf, style, order, M);

% Match MCMC exactly
sim.f = res.f;
sim.P = (res.P(1:N/2+1).*S(1:N/2+1))./max(res.P(1:N/2+1).*S(1:N/2+1));

% Convert to time domain
sim.pFreq = [sim.P conj(sim.P(end-1:-1:2))];
sim.pFreq(N/2+1) = real(sim.pFreq(N/2+1));
sim.pTime = ifft(sim.pFreq,'symmetric')/dt;

% Apply filtering
sim.pTimeFilt = bandpass_butterworth(sim.pTime, filterband, Fs, filterorder);

% Convert back to frequency domain
L = length(sim.pTimeFilt);
NFFT = 2^nextpow2(L);
sim.pFreqFilt = fft(sim.pTimeFilt, NFFT) / L;
f = Fs/2 * linspace(0,1,NFFT/2+1);

simSpec = abs(sim.pFreqFilt(1:NFFT/2+1));
simSpec = simSpec / max(simSpec);

%  8: Return PHYSICAL parameters, not normalized
globalBestIndiv = params_final;

end


