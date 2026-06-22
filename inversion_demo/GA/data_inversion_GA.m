%% INVERSION ETNA - GENETIC ALGORITHM VERSION %%
%
% Script file that performs GA inversion to invert harmonic infrasound
% observations for crater geometry using real data from Mount Etna
%
% Modified from MCMC version by Leighton Watson
% GA implementation by Daniel Spencer
% dsp65@uclive.ac.nz // daniel.spencer2007@gmail.com

%% Self-contained setup
scriptDir = fileparts(mfilename('fullpath'));
cd(scriptDir);
projectRoot = fullfile(scriptDir, '..');
restoredefaultpath;
rehash toolboxcache;
addpath(genpath(projectRoot));

%% Housekeeping
clearvars;
clc;
cmap = get(gca,'ColorOrder');
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
addpath(fullfile(rootDir, 'source', 'resonance'));
addpath(fullfile(rootDir, 'source', 'SBPoperators'));
addpath(fullfile(rootDir, 'source', 'inv'));

% Clear figures 1-3
for i = 1:3
    figure(i); clf;
end

%% Logicals %%
save_output = 0; % logical that determines if outputs are saved or not
plot_output = true; % logical that determines if outputs are plotted

%% Discretization parameters 
T = 25; 
N = 250; 
dt = T/N; 
Nyquist = 1/(2*dt); 
Nf = N/2+1;
discrParams = [T N Nf Nyquist dt];

%% Temperature parameters 
craterTemp = 100; 
atmoTemp = 0; 
temp = [craterTemp, atmoTemp]; 
M = problemParametersInv(craterTemp, atmoTemp);

%% Load real data
dataLoc = 'GRL2020\data';
dataStr = 'Etna2018Phase3 ';
S = load(fullfile(rootDir, dataLoc, dataStr));
if isfield(S,'dataF') && isfield(S,'dataAmp')
    dataF = S.dataF;
    dataAmp = S.dataAmp;
else
    error('MAT file does not contain dataF and dataAmp.');
end
data_freq = [dataF, dataAmp]; % format data

%% Filter properties
filterband = [0.25 4.8];
filterorder = 4; 
Fs = 10;
filterProps = [filterband, filterorder, Fs];
freqLim = 3;

%% GA Parameters
useTimeLimit = true;    % Boolean for constraint
maxTime = 60 * 30;      % seconds
popSize = 20;           % population size
nGen = 1800;            % number of generations (only used if useTimeLimit=false)
mutationRate = 0.14;    % probability of mutation
crossoverRate = 0.7;    % probability of crossover
eliteCount = 1;         % number of elite individuals to preserve
tournamentSize = 2;     % tournament selection size

%% Geometry parameters 
geomFlag = 1; 
geomR0 = 80; 
geomDepth = 200; 
geomParams = [geomDepth geomR0 geomR0 geomR0 geomR0 geomR0];
geomLowerBnds = [50 80 1 1 1 1];
geomUpperBnds = [300 150 140 120 120 120];
nx = length(geomParams)-1;

%% Source parameters
srcFlag = 0; 
srcParams = 0.3; 
srcStyle = 'Brune';
srcUpperBnds = 5; 
srcLowerBnds = 0.01;

%% Format parameters
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

%% Run GA inversion
tic
[population, misfit, simSpec, f, bestIndiv, globalBestIndiv, ...
        finalGen] = ga_spec(popSize, nGen, mutationRate, crossoverRate, ...
        eliteCount, tournamentSize, geomParams, geomFlag, srcParams, ...
        srcFlag, srcStyle, lowerBnds, upperBnds, discrParams, temp, ...
            filterProps, data_freq, freqLim, useTimeLimit, maxTime);
elapsed_time = toc;

fprintf('\nGA Inversion Complete:\n');
fprintf('  Total time: %.2f minutes\n', elapsed_time/60);
fprintf('  Generations completed: %d\n', finalGen);
fprintf('  Final best misfit: %.6f\n', min(misfit));

%% Compute best geometry
bestParams = lowerBnds + globalBestIndiv .* (upperBnds - lowerBnds);
shapeF = geomFunction(bestParams);
depthF = shapeF(1,1);

% Also compute initial geometry for comparison
shapeI = geomFunction(geomParams);
depthI = shapeI(1,1);
%% Plot outputs
if plot_output

    col.init  = [0 0.4470 0.7410];   % blue for initial geometry
    col.data  = col.init;            % blue for data spectrum
    col.mean  = [0.9290 0.6940 0.1250]; % yellow for post burn-in mean
    col.best  = [0.8500 0.3250 0.0980]; % red for best-misfit / best fit

    % Crater geometry comparison
    figure(1); clf; hold on;

    % Initial geometry
    plot(shapeI(:,2), shapeI(:,1), 'Color', col.init, 'LineWidth', 2);
    plot(-shapeI(:,2), shapeI(:,1), 'Color', col.init, 'LineWidth', 2);
    plot([-shapeI(1,2) shapeI(1,2)], [depthI depthI], 'Color', col.init, 'LineWidth', 2);

    % Post burn-in mean geometry
    burn_in = max(2, ceil(finalGen / 10));
    pop_trunc = population(burn_in:end, :);   % truncated population
    pop_physical = lowerBnds(1:length(geomParams)) + pop_trunc(:, 1:length(geomParams)) ...
        .* (upperBnds(1:length(geomParams)) - lowerBnds(1:length(geomParams)));
    valid = all(isfinite(population), 2);
    geomParams_mean = mean(population(valid, :), 1);
    shapeM = geomFunction(geomParams_mean);
    depthM = shapeM(1,1);
    
    plot(shapeM(:,2), shapeM(:,1), 'Color', col.mean, 'LineWidth', 2);
    plot(-shapeM(:,2), shapeM(:,1), 'Color', col.mean, 'LineWidth', 2);
    plot([-shapeM(1,2) shapeM(1,2)], [depthM depthM], 'Color', col.mean, 'LineWidth', 2);

    % Best-fit geometry
    shapeB = geomFunction(bestParams);
    depthB = shapeB(1,1);

    plot(shapeB(:,2), shapeB(:,1), 'Color', col.best, 'LineWidth', 2);
    plot(-shapeB(:,2), shapeB(:,1), 'Color', col.best, 'LineWidth', 2);
    plot([-shapeB(1,2) shapeB(1,2)], [depthB depthB], 'Color', col.best, 'LineWidth', 2);

    set(gca,'YDir','Reverse'); 
    xlabel('Radius (m)'); 
    ylabel('Depth (m)');
    title('Crater Geometry: GA Inversion');
    legend('Initial','','','Post burn-in mean','','','Best Fit','','','Location','best');
    axis equal; 
    grid on;

    % Spectra comparison
    figure(2); clf; hold on;

    % Data spectrum
    plot(dataF, abs(dataAmp)./max(abs(dataAmp)), 'Color', col.data, 'LineWidth', 2);

    % Post burn-in mean spectrum
    simSpec_trunc = simSpec(burn_in:end,:);
    spec_mean = mean(abs(simSpec_trunc),1);
    spec_mean = spec_mean ./ max(spec_mean);
    plot(f, spec_mean, 'Color', col.mean, 'LineWidth', 2);

    % Best-fit spectrum
    spec_best = abs(simSpec(end,:));
    spec_best = spec_best ./ max(spec_best);
    plot(f, spec_best, 'Color', col.best, 'LineWidth', 2);

    xlim([0 3]);
    xlabel('Frequency (Hz)'); 
    ylabel('Normalized Amplitude Spectra');
    title('Data vs GA Spectra');
    legend('Data','Post burn-in mean','Best Fit','Location','best'); 
    grid on;

    % Best-so-far misfit evolution
    figure(3); clf; hold on;
    bestMisfit = cummin(misfit);
    plot(bestMisfit, 'Color', col.best, 'LineWidth', 2);
    xline(burn_in, 'r--', 'LineWidth', 2);
    xlabel('Generation');
    ylabel('Best-so-far Misfit (L2 norm)');
    title('GA Misfit Evolution');
    xlim([2 finalGen]);
    grid on;

    % Print best parameters
    fprintf('\nBest fit parameters:\n');
    fprintf('  Depth: %.2f m\n', bestParams(1));
    for i = 2:length(bestParams)
        fprintf('  Radius %d: %.2f m\n', i-1, bestParams(i));
    end

end

%% Save outputs
if save_output
    pathname = fullfile(homeDir, 'OneDrive', 'Summer Research', 'Saved Results GA', 'Data Inversion Results');
    if ~exist(pathname, 'dir')
        mkdir(pathname);
    end
    filename = ['GA_InvOut_', dataStr, ...
                '_Ngen', num2str(finalGen), ...
                '_R0', num2str(geomR0), ...
                '_D', num2str(geomDepth), ...
                '_T', num2str(craterTemp), 'C', ...
                '_freqLim', num2str(freqLim), ...
                '_nx', num2str(nx), '.mat'];
    save(fullfile(pathname,filename), 'population', 'misfit', 'simSpec', 'f', ...
         'bestIndiv', 'globalBestIndiv', 'bestParams', 'finalGen', ...
         'params', 'geomParams', 'srcParams', 'geomFlag', 'srcFlag', 'srcStyle', ...
         'lowerBnds', 'upperBnds', 'freqLim', 'discrParams', 'M', 'filterProps', ...
         'elapsed_time', 'popSize', 'mutationRate', 'crossoverRate');
    fprintf('Results saved to: %s\n', fullfile(pathname,filename));
end
bestParams
