%% INVERSION ETNA %%
%
% Script file that performs MCMC inversion to invert harmonic infrasound
% observations for crater geometry using real data from Mount Etna
%
% Written by Leighton Watson
% March 10, 2020
%
% Time limiter added by Daniel Spencer
% January
% leightonwatson@stanford.edu // leightonmwatson@gmail.com


%% Self-contained setup (from original working script)
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
rootDir = fullfile(homeDir, 'OneDrive', 'Summer Research', 'Github Files');

% Add all required folders
addpath(fullfile(rootDir, 'GRL2020', 'data'));
addpath(fullfile(rootDir, 'inversion_demo'));
addpath(fullfile(rootDir, 'inversion_demo', 'GA'));
addpath(fullfile(rootDir, 'source', 'resonance'));
addpath(fullfile(rootDir, 'source', 'SBPoperators'));
addpath(fullfile(rootDir, 'source', 'inv'));

% Clear figures 1-3
for i = 1:5
    figure(i); clf;
end

save_output = 0; % logical that determines if outputs are saved or not
plot_output = 1; % logical that determines if outputs are plotted

%% resonance 1D %%

% set the parameters for the resonance1d calculations 
T = 25;                             % total time (s)
N = 250;                            % number of grid points (formulas assume even N)
dt = T/N;                           % time step (s)
Nyquist = 1/(2*dt);                 % Nyquist frequency (Hz)
Nf = N/2+1;                         % number of frequency samples
freq = [0 Nyquist];                 % frequency range (Hz)
discrParams = [T N Nf Nyquist dt];  % save parameters into array

craterTemp = 100;   % crater temperature
atmoTemp = 0;       % atmospheric temperature
temp = [craterTemp, atmoTemp]; 

order = 4;                                      % order of numerical scheme (4, 6 or 8)
style = 'baffled piston';                       % acoustic radiation model ('monopole' or ' baffled piston')
M = problemParametersInv(craterTemp,atmoTemp);  % problem parameters required for resonance1d

%% data %%

dataStr = 'Etna2018Phase1';
datafile = strcat(dataStr,'.mat');
load(datafile);                 % load data
data_freq = [dataF, dataAmp];   % format data

filterband = [0.25 4.8];        % frequency band to filter
filterorder = 4;                % order of butterworth filter
Fs = 10;                        % sampling frequency
filterProps = [filterband, filterorder, Fs]; % filter properties - same as for data

%% inversion parameters %%

useTimeLimit = true;    % Boolean for time constraint
maxTime = 60 * 30;        % maximum time in seconds

nIter = 10000; % number of steps (only used if useTimeLimit=false)

dx = 0.1; % step size % use step size of 0.05 for paper inversions

freqLim = 3; % high cut frequency limit for misfit function (Hz)

%%% geometry parameters %%%
geomFlag = 1; % invert for geometry (boolean, 0 = no, 1 = yes)
geomR0 = 100; % radius of initial cylinder
geomDepth = 200; % depth 
geomParams = [geomDepth geomR0 geomR0 geomR0 geomR0 geomR0]; % first value is depth, other values are radius points that are equally spaced
geomLowerBnds = [50 80 1 1 1 1];
geomUpperBnds = [300 150 140 120 120 120];
nx = length(geomParams)-1; % number of geometry parameters

%%% source parameters %%%
srcFlag = 0; %invert for source (boolean, 0 = no, 1 = yes)
srcParams = 0.3;
srcStyle = 'Brune';
srcUpperBnds = [5];
srcLowerBnds = [0.01];

%%% format parameters %%%
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


%% mcmc inversion %%

tic % start timing
[x, misfit, simSpec, f, count] = mcmc_spec(nIter, dx, ... % perform MCMC inversion
    geomParams, geomFlag, srcParams, srcFlag, srcStyle,...
    lowerBnds, upperBnds, discrParams, temp, ...
    filterProps, data_freq, freqLim, useTimeLimit, maxTime);
elapsed_time = toc; % finish timing

disp(['Elapsed time is ',num2str(elapsed_time/60), ' minutes.']); % display timing in minutes

%% save outputs %%

% format spectra
spec.int = simSpec(1,:); % initial spectra
spec.fin = simSpec(end,:); % final spectra
burn_in = ceil(count/10); % remove the first 10% of successful samples to remove effect of initial conditions
simSpec_trunc = simSpec(burn_in:end,:); % remove burn in to reduce sensitivity to initial conditions
simSpec_mean = mean(abs(simSpec_trunc)); % mean spectra
spec.mean = abs(simSpec_mean)./max(abs(simSpec_mean)); % normalize spectra

pathname = strcat(pwd,'/invOutput'); % directory to save outputs
if useTimeLimit
    filename = strcat('DataInvOut','_',dataStr,'_',... % file name
        'Time',num2str(maxTime/60),'min_',...
        'count',num2str(count),'_',...
        'R0',num2str(geomR0),'m_','D',num2str(geomDepth),'m_',...
        'T',num2str(craterTemp),'C_',...
        'freqLim',num2str(freqLim),'_nx',num2str(nx));
else
    filename = strcat('DataInvOut','_',dataStr,'_',... % file name
        'Nit',num2str(nIter),'_',...
        'R0',num2str(geomR0),'m_','D',num2str(geomDepth),'m_',...
        'T',num2str(craterTemp),'C_',...
        'freqLim',num2str(freqLim),'_nx',num2str(nx));
end
matfile = fullfile(pathname,filename); % path and file for saving outputs

if save_output == 1
    save(matfile,'x','misfit','spec','f','count',...
    'params','geomParams','srcParams','geomFlag','srcFlag','srcStyle',...
    'lowerBnds','upperBnds','nIter','dx','freqLim',...
    'discrParams','M','filterProps','elapsed_time','useTimeLimit','maxTime');
end

%% plot outputs %%

if plot_output == 1
    
    burn_in = ceil(count/10);

    % Post burn-in samples
    x_trunc = x(burn_in:end,:);
    simSpec_trunc = simSpec(burn_in:end,:);
    
    % Best-misfit index
    [~, bestIdx] = min(misfit);
    
    col.init = [0 0.4470 0.7410];     % blue for initial geometry
    col.data = col.init;              % blue for data
    col.mean = [0.9290 0.6940 0.1250]; % yellow
    col.best = [0.8500 0.3250 0.0980]; % red
    figure(1); clf; hold on
    
    %%% Initial estimate of crater geometry
    shapeI = geomFunction(geomParams);
    depthI = shapeI(1,1);
    
    plot( shapeI(:,2), shapeI(:,1), 'Color', col.init );
    plot(-shapeI(:,2), shapeI(:,1), 'Color', col.init );
    plot([-shapeI(1,2) shapeI(1,2)], [depthI depthI], 'Color', col.init );
    
    %%% Post burn-in mean crater geometry
    x_mean = mean(x_trunc(:,1:length(geomParams)));
    geomParams_mean = x_mean .* geomParams;
    shapeM = geomFunction(geomParams_mean);
    depthM = shapeM(1,1);
    plot( shapeM(:,2), shapeM(:,1), 'Color', col.mean );
    plot(-shapeM(:,2), shapeM(:,1), 'Color', col.mean );
    plot([-shapeM(1,2) shapeM(1,2)], [depthM depthM], 'Color', col.mean );
    
    %%% Best misfit geometry %%%
    x_best = x(bestIdx,1:length(geomParams));
    geomParams_best = x_best .* geomParams;

    shapeB = geomFunction(geomParams_best);
    depthB = shapeB(1,1);

    plot( shapeB(:,2), shapeB(:,1), 'Color', col.best );
    plot(-shapeB(:,2), shapeB(:,1), 'Color', col.best );
    plot([-shapeB(1,2) shapeB(1,2)], [depthB depthB], 'Color', col.best );

    set(gca,'YDir','Reverse');
    xlabel('Radius (m)');
    ylabel('Depth (m)');
    title('Crater Geometry: MCMC Inversion');

    legend('Initial','','', ...
           'Post burn-in mean','','', ...
           'Best misfit','','', ...
           'Location','best');
    grid on;
    
    %%% spectra %%%

    figure(2); clf; hold on;
    
    %%% Data
    plot(dataF, abs(dataAmp)./max(abs(dataAmp)), ...
         'Color', col.data, 'LineWidth', 2);
    
    
    %%% Post burn-in mean estimate of spectra
    spec_mean = mean(abs(simSpec_trunc),1);
    spec_mean = spec_mean ./ max(spec_mean);

    plot(f, spec_mean, ...
         'Color', col.mean, 'LineWidth', 2);
    
    %%% Best misfit spectra
    
    spec_best = abs(simSpec(bestIdx,:));
    spec_best = spec_best ./ max(spec_best);

    plot(f, spec_best, ...
         'Color', col.best, 'LineWidth', 2);

    xlim([0 3]);
    xlabel('Frequency (Hz)');
    ylabel('Normalized Amplitude Spectra');
    title('Observed vs MCMC Spectra');

    legend('Data', 'Post burn-in mean', 'Best misfit', 'Location','best');
    grid on;

    %%% Best current misfit evolution %%%
    figure(3); clf;
    bestMisfit = cummin(misfit);
    
    plot(bestMisfit, 'LineWidth', 2); hold on;
    xline(burn_in, 'r--', 'LineWidth', 2);
    
    xlabel('Successful Iteration');
    ylabel('Best Current Misfit (L2 norm)');
    title('Best Misfit Evolution');
    
    xlim([1 count]);
    grid on

    %%% parameter histograms (post burn-in) %%%
    figure(4); clf;
    k = length(params);
    for i = 1:k
        paramPlot = params(i) * x_trunc(:,i);

        subplot(1,k,i); hold on; box on;
        histogram(paramPlot);

        xline(lowerBnds(i), 'r--', 'LineWidth', 2);
        xline(upperBnds(i), 'r--', 'LineWidth', 2);

        xlim([lowerBnds(i) upperBnds(i)]);
        xlabel(['Param ' num2str(i)]);
    end
    sgtitle('Parameter Distributions (post burn-in)');

	%%% parameter time evolution %%%
    figure(5); clf;

    for i = 1:k
        paramPlot = params(i) * x_trunc(:,i);

        subplot(1,k,i); hold on; box on;
        plot(paramPlot);
        xline(burn_in, 'r--', 'LineWidth', 2);

        xlabel('Iteration');
        ylabel(['Param ' num2str(i)]);
    end
    sgtitle('Parameter Evolution (post burn-in)');
    
    % Print summary statistics
    fprintf('\n=== MCMC Summary Statistics ===\n');
    fprintf('Mean parameter values (post burn-in):\n');
    stDev = [0, 0, 0, 0, 0, 0];
    for i = 1:k
        paramPlot = params(i)*x_trunc(:,i);
        stDev(i) = std(paramPlot);
        fprintf('  Param %d: %.2f ± %.2f\n', i, mean(paramPlot), std(paramPlot));
    end
    
end
min(misfit)
mean(misfit(burn_in:end))
geomParams_best
geomParams_mean
stDev