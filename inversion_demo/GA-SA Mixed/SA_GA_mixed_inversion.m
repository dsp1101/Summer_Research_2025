%% INVERSION ETNA - MIXED METHOD IMPLEMENTATION %%
%
% Script file that performs SA to identify minima, then GA to optimise
% within the minima, to invert harmonic infrasound observations for crater
% geometry using real data from Mount Etna
%
% Modified from SA and GA versions by Coleman Campbell and Daniel Spencer
% by Daniel Spencer
%
% dsp65@uclive.ac.nz // daniel.spencer2007@gmail.com
% 23/01/2026

clearvars, close all hidden; clc;
echo off;

set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',18);

% Reset path
path(pathdef)
base = userpath;

% Build user-independent root folder
homeDir = getenv('USERPROFILE');
rootDir = fullfile(homeDir, 'OneDrive', 'Summer Research 25-26', 'Github Files');

% Add all required folders
addpath(fullfile(rootDir, 'GRL2020', 'data'));
addpath(fullfile(rootDir, 'inversion_demo'));
addpath(fullfile(rootDir, 'inversion_demo', 'GA'));
addpath(fullfile(rootDir, 'inversion_demo', 'GA-SA Mixed'));
addpath(fullfile(rootDir, 'source', 'resonance'));
addpath(fullfile(rootDir, 'source', 'SBPoperators'));
addpath(fullfile(rootDir, 'source', 'inv'));

save_output = 1;
plot_output = 1;

%% SA inversion parameters
max_time_minutes = 10;
max_time = max_time_minutes * 60; 
freqLim = 3;

%% CRATER SETUP
craterTemp = 100;
atmoTemp = 0;
temp = [craterTemp, atmoTemp];

% geometry parameters
geomFlag = 1;
geomR0 = 100;
geomDepth = 200;
geomParams = [geomDepth geomR0 geomR0 geomR0 geomR0 geomR0];
geomLowerBnds = [50 80 1 1 1 1];
geomUpperBnds = [300 150 140 120 120 120];
nx = length(geomParams)-1;

%% LOAD DATA
dataStr = 'Etna2018phase1';
datafile = strcat(dataStr,'.mat');
load(datafile);
data_freq = [dataF, dataAmp];

%% miscellaneous setup
T = 25;
N = 250;
dt = T/N;
Nyquist = 1/(2*dt);
Nf = N/2+1;
freq = [0 Nyquist];
discrParams = [T N Nf Nyquist dt];

order = 4;
style = 'baffled piston';
M = problemParametersInv(craterTemp,atmoTemp);

filterband = [0.25 4.8];
filterorder = 4;
Fs = 10;
filterProps = [filterband, filterorder, Fs];

% source parameters
srcFlag = 0;
srcParams = 0.3;
srcStyle = 'Brune';
srcUpperBnds = 5;
srcLowerBnds = 0.01;

% format parameters
if geomFlag && srcFlag
    upperBnds = [geomUpperBnds srcUpperBnds];
    lowerBnds = [geomLowerBnds srcLowerBnds];
    params = [geomParams srcParams];
elseif geomFlag
    upperBnds = geomUpperBnds;
    lowerBnds = geomLowerBnds;
    params = geomParams;
elseif srcFlag
    upperBnds = srcUpperBnds;
    lowerBnds = srcLowerBnds;
    params = srcParams;
end

%% create handle class for tracking history
estimatedIterations=ceil(max_time*12);
history = OptimisationHistory(estimatedIterations);

%% Get objective function for simulannealbnd
objectiveFcn = @(x_norm) SA_compute_misfit(x_norm, params, geomParams, ... 
    geomFlag, srcParams, srcFlag, srcStyle, discrParams, temp, ...
    filterProps, data_freq, freqLim);

%% configure Simulated Annealing options
outputFcn = @(options, optimvalues, flag) trackProgress(options, ...
    optimvalues, flag, history);

options = optimoptions('simulannealbnd', ...
    'AnnealingFcn', @annealingboltz, ...
    'MaxTime', max_time, ...
    'maxIterations', Inf, ...
    'OutputFcn', outputFcn, ...
    'Display', 'final', ...
    'InitialTemperature', 150, ...
    'TemperatureFcn', @temperatureexp, ...
    'ReannealInterval', 80, ...
    'MaxFunctionEvaluations', Inf, ...
    'FunctionTolerance', 0, ...
    'MaxStallIterations', 2147483647);

%% Initial guess
x0=ones(size(params));

% normalise bounds
lb_normalised = lowerBnds ./ params;
ub_normalised = upperBnds ./ params;

%% Run simulated annealing
disp('========================================');
disp('PHASE 1: SIMULATED ANNEALING');
disp('========================================');
tic;
[x_optimal, final_misfit, exitflag, output] = simulannealbnd(...
    objectiveFcn, x0, lb_normalised, ub_normalised, options);
elapsed_time=toc;

fprintf('SA completed in %.6f minutes\n', elapsed_time/60);
fprintf('Total iterations:  %d\n', output.iterations);
fprintf('Final misfit: %.6f\n', final_misfit);

%% Extract results from history
misfit = history.getMisfits();
x_normalized = history.getParams();  % These are normalized (relative to x0=ones)

% Convert normalized parameters to physical parameters
% SA explores around x0 = ones(size(params)), so:
% physical_params = x_normalized .* params
x_physical = x_normalized .* params;

params_optimal = x_optimal .* params;

%% Sort top local minima/basins
minWindowSize = 100;  
localMinIdx = islocalmin(misfit, 'MinSeparation', minWindowSize); 

% check if last basin is included
if length(misfit) > minWindowSize
    final_window_start = length(misfit) - minWindowSize + 1;
    final_window = misfit(final_window_start:end);
    [final_min_val, final_min_rel_idx] = min(final_window);
    final_min_abs_idx = final_window_start + final_min_rel_idx - 1;
    
    if ~localMinIdx(final_min_abs_idx)
        localMinIdx(final_min_abs_idx) = true;
    end
end

localMinValues = misfit(localMinIdx);
localMinIndices = find(localMinIdx);
localMinParams = x_physical(localMinIndices, :); 

% Sort local minima by value
[sortedMinVals, sortOrder] = sort(localMinValues);
sortedMinIndices = localMinIndices(sortOrder);
sortedMinParams = localMinParams(sortOrder, :);

%% Initialize arrays with correct dimensions and consistent naming
nBestParams = 2;
best_misfits = zeros(nBestParams, 1);  % Column vector
best_indices = zeros(nBestParams, 1);  % Column vector
best_params = zeros(nBestParams, size(sortedMinParams, 2));  % Matrix

% Remove duplicates based on parameter similarity
nBasins = 0;
param_tolerance = 0.1;

for i = 1:length(sortedMinVals)
    if nBasins == 0
        % First basin
        nBasins = 1;
        best_misfits(1) = sortedMinVals(i);
        best_indices(1) = sortedMinIndices(i);
        best_params(1, :) = sortedMinParams(i, :);
    else
        % Check if this minimum is in a new basin
        is_new_basin = true;
        for j = 1:nBasins
            param_diff = abs(sortedMinParams(i, :) - best_params(j, :));
            if all(param_diff < param_tolerance)
                is_new_basin = false;
                break;
            end
        end
        
        if is_new_basin
            nBasins = nBasins + 1;
            best_misfits(nBasins) = sortedMinVals(i);
            best_indices(nBasins) = sortedMinIndices(i);
            best_params(nBasins, :) = sortedMinParams(i, :);
        end
    end
    
    if nBasins >= nBestParams
        break;
    end
end

% Trim to actual number found
best_misfits = best_misfits(1:nBasins);
best_indices = best_indices(1:nBasins);
best_params = best_params(1:nBasins, :);

%% FIX: Denormalize SA parameters before passing to GA
% Parameters are now already in physical units from x_denorm
% But we should verify the bounds make sense
best_params_physical = best_params;  % Already physical from x_denorm

% Verify parameters are within bounds
for i = 1:nBasins
    if any(best_params_physical(i,:) < lowerBnds) || any(best_params_physical(i,:) > upperBnds)
        warning('Basin %d parameters outside bounds. Clipping...', i);
        best_params_physical(i,:) = max(lowerBnds, min(upperBnds, best_params_physical(i,:)));
    end
end

disp('========================================');
disp('PHASE 2: GENETIC ALGORITHM REFINEMENT');
disp('========================================');
fprintf('Found %d distinct basins from SA\n', nBasins);
fprintf('Running GA on each basin...\n');

%% GA Parameters
useTimeLimit = true;
maxTime = 60 * 10;  % seconds per basin
popSize = 24;
nGen = 1800;
mutationRate = 0.18;
crossoverRate = 0.7;
eliteCount = 1;
tournamentSize = 3;
localSearchRadius = [5, 2.5, 2.5, 2.5, 2.5, 2.5];  % Radius for [depth, r1, r2, r3, r4, r5]

% Output data holders
finalParams = zeros(nBasins, size(best_params_physical, 2));
finalSpecs = [];  % Initialize empty, will grow dynamically
finalMisfits = zeros(nBasins, 1);
f_common = [];  % Store frequency vector from first run

for i = 1:nBasins
    fprintf('\n--- Refining Basin %d/%d (SA misfit: %.6f) ---\n', ...
            i, nBasins, best_misfits(i));
    
    initialParams = best_params_physical(i, :);  % Use physical parameters
    lowerBnds_local = initialParams - localSearchRadius;
    upperBnds_local = initialParams + localSearchRadius;
    
    % Clip to global bounds (can't go outside SA's explored region)
    lowerBnds_local = max(lowerBnds_local, lowerBnds);
    upperBnds_local = min(upperBnds_local, upperBnds);
    
    % Ensure lower < upper (in case initial param is at boundary)
    for j = 1:length(lowerBnds_local)
        if lowerBnds_local(j) >= upperBnds_local(j)
            % Give some wiggle room
            mid = (lowerBnds(j) + upperBnds(j)) / 2;
            lowerBnds_local(j) = max(lowerBnds(j), mid - 5);
            upperBnds_local(j) = min(upperBnds(j), mid + 5);
        end
    end
    
    fprintf('Local bounds for this basin:\n');
    fprintf('  Lower: ');
    fprintf('%.2f ', lowerBnds_local);
    fprintf('\n  Upper: ');
    fprintf('%.2f ', upperBnds_local);
    fprintf('\n');
    %% Run GA inversion
    tic
    [population, misfit_ga, simSpec, f, bestIndiv, globalBestIndiv, ...
        finalGen] = ga_local(popSize, nGen, mutationRate, crossoverRate, ...
        eliteCount, tournamentSize, geomParams, geomFlag, srcParams, ...
        srcFlag, srcStyle, lowerBnds, upperBnds, discrParams, temp, ...
        filterProps, data_freq, freqLim, useTimeLimit, maxTime, ...
        initialParams);
    elapsed_time = toc;
    
    % Store results
    finalParams(i, :) = globalBestIndiv;
    
    % Handle variable spectrum length
    if i == 1
        f_common = f;
        finalSpecs = zeros(nBasins, length(simSpec));
    end
    
    % Interpolate if frequency vectors don't match
    if length(f) ~= length(f_common)
        warning('Basin %d has different frequency vector length. Interpolating...', i);
        simSpec = interp1(f, simSpec, f_common, 'pchip', 0);
    end
    
    finalSpecs(i, :) = simSpec;
    finalMisfits(i) = misfit_ga(end);
    
    fprintf('Basin %d completed in %.2f min, final misfit: %.6f\n', ...
            i, elapsed_time/60, finalMisfits(i));
end

%% Find overall best solution
[overallBestMisfit, bestBasinIdx] = min(finalMisfits);
fprintf('\n========================================\n');
fprintf('OPTIMIZATION COMPLETE\n');
fprintf('========================================\n');
fprintf('Best solution from Basin %d\n', bestBasinIdx);
fprintf('Final misfit: %.6f\n', overallBestMisfit);
fprintf('Best parameters: ');
fprintf('%.2f ', finalParams(bestBasinIdx, :));
fprintf('\n');

%% Plotting
if plot_output
    colors = lines(nBasins + 1);
    col_init = colors(1,:);
    col_basins = colors(2:end,:);
    
    % Plot Geometries
    figure(1); clf; hold on;
    shapeI = geomFunction(geomParams);
    depthI = shapeI(1,1);
    plot(shapeI(:,2), shapeI(:,1), 'Color', col_init, 'LineWidth', 2);
    plot(-shapeI(:,2), shapeI(:,1), 'Color', col_init, 'LineWidth', 2);
    plot([-shapeI(1,2) shapeI(1,2)], [depthI depthI], 'Color', col_init, 'LineWidth', 2);
    
    for i = 1:nBasins
        paramsA = finalParams(i, :);  % Use GA-refined parameters
        shapeA = geomFunction(paramsA);
        depthA = shapeA(1,1);
        col = col_basins(i,:);
        
        % Add transparency for non-best solutions
        if i ~= bestBasinIdx
            col = [col 0.4];
        end
        
        plot(shapeA(:,2), shapeA(:,1), 'Color', col, 'LineWidth', 2);
        plot(-shapeA(:,2), shapeA(:,1), 'Color', col, 'LineWidth', 2);
        plot([-shapeA(1,2) shapeA(1,2)], [depthA depthA], 'Color', col, 'LineWidth', 2);
    end
    set(gca,'YDir','Reverse'); 
    xlabel('Radius (m)'); 
    ylabel('Depth (m)');
    title('Crater Geometry: SA-GA Hybrid Inversion');
    legend('Initial', '', '', 'Basin solutions', 'Location','best');
    axis equal; 
    grid on;
    
    % Plot Spectra
    figure(2); clf; hold on;
    plot(dataF, abs(dataAmp)./max(abs(dataAmp)), 'k', 'LineWidth', 3);
    spec = abs(finalSpecs(bestBasinIdx,:));
    spec = spec ./ max(spec);
    plot(f_common, spec, 'Color', col_basins(bestBasinIdx,:), 'LineWidth', 2.5);
    
    xlim([0 3]);
    xlabel('Frequency (Hz)');
    ylabel('Normalized Amplitude Spectrum');
    title('Data vs Simulated Spectra (SA-GA Hybrid)');
    legend({'Data', 'Best fit'}, 'Location','best');
    grid on;
    
    % Plot SA misfit evolution
    figure(3); clf;
    plot(misfit, 'LineWidth', 1.5);
    hold on;
    plot(best_indices, best_misfits, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('SA Iteration');
    ylabel('Misfit');
    title('SA Exploration with Identified Basins');
    legend('Misfit trajectory', 'Selected basins', 'Location', 'best');
    grid on;
end