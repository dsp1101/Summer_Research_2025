%% SIMULATED ANNEALING INVERSION %%
% Script file that performs simmulated annealing inversion to invert
% harmonic infrasound observations for crater geometry using real data from 
% Mount Etna. Simulated inversion 
% 
% Set up to run for max_time_minutes
% 
% Applicable to other craters if setup geometries and conditions altered, 
% as found in sections _CRATER SETUP_ and _LOAD DATA_
% 

clearvars, close all hidden; clc;
echo off;

cmap = get(gca,'ColorOrder');
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',18);

%***change paths depending on CRes folder structure***
path(pathdef);
addpath ..\CResPapers\GRL2020\data\
addpath ..\..\source\inv\
addpath ..\..\source\resonance\
addpath ..\..\source\SBPoperators\


save_output = 1; % logical that determines if outputs are saved or not
plot_output = 1; % logical that determines if outputs are plotted


%% inversion parameters
max_time_minutes = 600; % minutes to run code
max_time = max_time_minutes*60; % convert to seconds
freqLim = 3; % high cut frequency limit for misfit function (Hz)

%% CRATER SETUP %%
% Tweak parameters depending on target crater. Current set up for Mount
% Etna, in particular for 2018EtnaPhase1 and 2018EtnaPhase3 data (GRL2020)

craterTemp = 100; % crater temperature
atmoTemp = 0; % atmospheric temperature
temp = [craterTemp, atmoTemp]; % Celsius

% geometry parameters
geomFlag = 1; % invert for geometry (boolean, 0 = no, 1 = yes)
geomR0 = 100; % radius of initial cylinder. See geomLowerBnds and geomUpperBnds
geomDepth = 200; % depth. See geomLowerBnds and geomUpperBnds
geomParams = [geomDepth geomR0 geomR0 geomR0 geomR0 geomR0]; % first value is depth, other values are radius points that are equally spaced
geomLowerBnds = [50 80 1 1 1 1]; %[50 80 1 1 1 1] for Etna
geomUpperBnds = [300 150 140 120 120 120]; %[300 120 120 120 120 120] for Etna
nx = length(geomParams)-1; % number of geometry parameters


%% LOAD DATA (sample) %%
dataStr = 'Etna2018phase1';
datafile = strcat(dataStr,'.mat');
load(datafile); % load data
data_freq = [dataF, dataAmp]; % format data

%% miscellaneous setup %%

% resonance 1D - set the parameters for the resonance1d calculations 
T = 25; % total time (s)
N = 250; % number of grid points (formulas assume even N)
dt = T/N; % time step (s)
Nyquist = 1/(2*dt); % Nyquist frequency (Hz)
Nf = N/2+1; % number of frequency samples
freq = [0 Nyquist]; % frequency range (Hz)
discrParams = [T N Nf Nyquist dt]; % save parameters into array

order = 4; % order of numerical scheme (4, 6 or 8)
style = 'baffled piston'; % acoustic radiation model ('monopole' or ' baffled piston')
M = problemParametersInv(craterTemp,atmoTemp); % problem parameters required for resonance1d

filterband = [0.25 4.8]; % frequency band to filter
filterorder = 4; % order of butterworth filter
Fs = 10; % sampling frequency
filterProps = [filterband, filterorder, Fs]; % filter properties - same as for data

% source parameters %
srcFlag = 0; %invert for source (boolean, 0 = no, 1 = yes)
srcParams = 0.3;
srcStyle = 'Brune';
srcUpperBnds = 5;
srcLowerBnds = 0.01;

% format parameters %
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

%% create handle class for tracking history %%
estimatedIterations=ceil(max_time*12); % estimate no more than 12 iterations per second, average
history = OptimisationHistory(estimatedIterations);

%% Get objective function for simulannealbnd

objectiveFcn = @(x_norm) SA_compute_misfit(x_norm, params, geomParams, ... 
    geomFlag, srcParams, srcFlag, srcStyle, discrParams, temp, ...
    filterProps, data_freq, freqLim);

%% configure Simulated Annealing options %%
% create output function 
outputFcn = @(options, optimvalues, flag) trackProgress(options, ...
    optimvalues, flag, history);

options = optimoptions('simulannealbnd', ...
    'AnnealingFcn', @annealingboltz, ... % step length = sqrt of temperature
    'MaxTime', max_time, ...
    'maxIterations', Inf, ...
    'OutputFcn', outputFcn, ...
    'Display', 'final', ... % custom display in trackProgress
    'InitialTemperature', 100, ...
    'TemperatureFcn', @temperatureexp, ...
    'ReannealInterval', 100, ...
    'MaxFunctionEvaluations', Inf, ...
    'FunctionTolerance', 0, ...
    'MaxStallIterations', 2147483647);

%% Initial guess
x0=ones(size(params));

% normalise bounds
lb_normalised = lowerBnds ./ params;
ub_normalised = upperBnds ./ params;

%% Run simulated annealing
tic;
[x_optimal, final_misfit, exitflag, output] = simulannealbnd(...
    objectiveFcn, x0, lb_normalised, ub_normalised, options);
elapsed_time=toc;

disp(strcat('Optimisation completed in ', num2str(elapsed_time/60), ' minutes.'));
disp(strcat('Total iterations: ',num2str(output.iterations)));
disp(strcat('Final misfit: ',num2str(final_misfit)));

%% Extract results from history (handle class)
misfit = history.getMisfits();
x = history.getParams();

params_optimal = x_optimal.* params;

%% Sort top local minima/basins
minWindowSize = 100;  
localMinIdx = islocalmin(misfit, 'MinSeparation', minWindowSize); 

% check if last basin is included
if length(misfit) > minWindowSize
    % Look at last minWindowSize points
    final_window_start = length(misfit) - minWindowSize + 1;
    final_window = misfit(final_window_start:end);
    [final_min_val, final_min_rel_idx] = min(final_window);
    final_min_abs_idx = final_window_start + final_min_rel_idx - 1;
    
    % Check if this minimum was already detected
    if ~localMinIdx(final_min_abs_idx)
        % Add it as a local minimum
        localMinIdx(final_min_abs_idx) = true;
    end
end

localMinValues = misfit(localMinIdx);
localMinIndices = find(localMinIdx);
localMinParams = x(localMinIndices, :);

% Sort local minima by value
[sortedMinVals, sortOrder] = sort(localMinValues);
sortedMinIndices = localMinIndices(sortOrder);
sortedMinParams = localMinParams(sortOrder, :);

% Remove duplicates based on parameter similarity 
nBasins = 0;
best5_misfits = [];
best5_indices = [];
best5_params = [];

param_tolerance = 0.1; 
for i = 1:length(sortedMinVals)
    if nBasins == 0
        % First basin
        nBasins = 1;
        best5_misfits(nBasins) = sortedMinVals(i);
        best5_indices(nBasins) = sortedMinIndices(i);
        best5_params(nBasins, :) = sortedMinParams(i, :);
    else
        % Check if this minimum is in a new basin (different parameters)
        is_new_basin = true;
        for j = 1:nBasins
            param_diff = abs(sortedMinParams(i, :) - best5_params(j, :));
            if all(param_diff < param_tolerance)
                is_new_basin = false;
                break;
            end
        end
        
        if is_new_basin
            nBasins = nBasins + 1;
            best5_misfits(nBasins) = sortedMinVals(i);
            best5_indices(nBasins) = sortedMinIndices(i);
            best5_params(nBasins, :) = sortedMinParams(i, :);
        end
    end
    
    if nBasins >= 5
        break;
    end
end


fprintf('Best %d distinct basins:\n', nBasins);
for i = 1:nBasins
    fprintf('Basin %d: misfit = %.6f (local min at iteration %d), ', ...
        i, best5_misfits(i), best5_indices(i));
    fprintf('Parameters: [%s]\n', num2str(best5_params(i, :).*geomParams, '%.4f '));
end

%% Compute final spectrum for Optimal Parameters
% In order to plot 3 best basins
params_optimalA=best5_params(1,:).*geomParams;
params_optimalB=best5_params(2,:).*geomParams;
params_optimalC=best5_params(3,:).*geomParams;

if geomFlag
    shapes.shapeFA = geomFunction(params_optimalA(1:length(geomParams)));
    shapes.depthFA = shapes.shapeFA(1,1);
    shapes.shapeFB = geomFunction(params_optimalB(1:length(geomParams)));
    shapes.depthFB = shapes.shapeFB(1,1);
    shapes.shapeFC = geomFunction(params_optimalC(1:length(geomParams)));
    shapes.depthFC = shapes.shapeFC(1,1);
else
    shapes.shapeF = geomFunction(geomParams);
    shapes.depthF = shapes.shapeF(1,1);
end

if srcFlag
    if geomFlag
        [S, ~, ~, ~] = sourceFunction(1, params_optimal(length(geomParams)+1:end), srcStyle, discrParams);
    else
        [S, ~, ~, ~] = sourceFunction(1, params_optimal(1:length(srcParams)), srcStyle, discrParams);
    end
else
    [S, ~, ~, ~] = sourceFunction(1, srcParams, srcStyle, discrParams);
end

% data spectrum
dataSpectrum.f = data_freq(:,1);
dataSpectrum.P = abs(data_freq(:,2))./max(abs(data_freq(:,2)));
f = dataSpectrum.f;

% Crater B (second best)
res = resonance1d(shapes.shapeFB, shapes.depthFB, freq, Nf, style, order, M);
sim.f=res.f;
sim.P = (res.P(1:N/2+1).*S(1:N/2+1))./max(res.P(1:N/2+1).*S(1:N/2+1));
sim.pFreq = [sim.P conj(sim.P(end-1:-1:2))];
sim.pFreq(N/2+1) = real(sim.pFreq(N/2+1));
sim.pTime = ifft(sim.pFreq,'symmetric')/dt;
sim.pTimeFilt = bandpass_butterworth(sim.pTime,filterband,Fs,filterorder);
L = length(sim.pTimeFilt);
NFFT = 2^nextpow2(L);
sim.pFreqFilt = fft(sim.pTimeFilt,NFFT)/L;
sim.pFreqAbsNorm = abs(sim.pFreqFilt(1:NFFT/2+1))./max(abs(sim.pFreqFilt(1:NFFT/2+1)));
sim.fFilt = Fs/2*linspace(0,1,NFFT/2+1);
spec.finB = pchip(sim.fFilt,sim.pFreqAbsNorm,dataSpectrum.f)';

% Crater C (third best)
res = resonance1d(shapes.shapeFC, shapes.depthFC, freq, Nf, style, order, M);
sim.f=res.f;
sim.P = (res.P(1:N/2+1).*S(1:N/2+1))./max(res.P(1:N/2+1).*S(1:N/2+1));
sim.pFreq = [sim.P conj(sim.P(end-1:-1:2))];
sim.pFreq(N/2+1) = real(sim.pFreq(N/2+1));
sim.pTime = ifft(sim.pFreq,'symmetric')/dt;
sim.pTimeFilt = bandpass_butterworth(sim.pTime,filterband,Fs,filterorder);
L = length(sim.pTimeFilt);
NFFT = 2^nextpow2(L);
sim.pFreqFilt = fft(sim.pTimeFilt,NFFT)/L;
sim.pFreqAbsNorm = abs(sim.pFreqFilt(1:NFFT/2+1))./max(abs(sim.pFreqFilt(1:NFFT/2+1)));
sim.fFilt = Fs/2*linspace(0,1,NFFT/2+1);
spec.finC = pchip(sim.fFilt,sim.pFreqAbsNorm,dataSpectrum.f)';

% Crater A (optimal). Last so that sim contains values from optimal crater
res = resonance1d(shapes.shapeFA, shapes.depthFA, freq, Nf, style, order, M);
sim.f=res.f;
sim.P = (res.P(1:N/2+1).*S(1:N/2+1))./max(res.P(1:N/2+1).*S(1:N/2+1));
sim.pFreq = [sim.P conj(sim.P(end-1:-1:2))];
sim.pFreq(N/2+1) = real(sim.pFreq(N/2+1));
sim.pTime = ifft(sim.pFreq,'symmetric')/dt;
sim.pTimeFilt = bandpass_butterworth(sim.pTime,filterband,Fs,filterorder);
L = length(sim.pTimeFilt);
NFFT = 2^nextpow2(L);
sim.pFreqFilt = fft(sim.pTimeFilt,NFFT)/L;
sim.pFreqAbsNorm = abs(sim.pFreqFilt(1:NFFT/2+1))./max(abs(sim.pFreqFilt(1:NFFT/2+1)));
sim.fFilt = Fs/2*linspace(0,1,NFFT/2+1);
spec.finA = pchip(sim.fFilt,sim.pFreqAbsNorm,dataSpectrum.f)';



%% Save Outputs
datestr=string(datetime('now','Format','yyyy-MM-dd_HH_mm_ss'));
pathname = strcat(pwd,'/simulated_annealing_output_data');
filename = strcat('SAInvOut','_',dataStr,'_',...
    'runlength',num2str(max_time_minutes),'min_',...
    'date', datestr, '_', ...
    'R0',num2str(geomR0),'m_','D',num2str(geomDepth),'m_',...
    'T',num2str(craterTemp),'C_',...
    'freqLim',num2str(freqLim),'_nx',num2str(nx));
matfile = fullfile(pathname,filename);

if save_output == 1
    save(matfile,'x','misfit','spec','f','params_optimal',...
        'params','geomParams','srcParams','geomFlag','srcFlag','srcStyle',...
        'lowerBnds','upperBnds','max_time','freqLim',...
        'discrParams','M','filterProps','elapsed_time','output', ...
        'best5_misfits', 'best5_indices', 'best5_params', 'options',...
        'data_freq', 'shapes');
end

%% Plot Outputs
if plot_output == 1
    
    % Crater geometry
    figure(1);
    shapeI = geomFunction(geomParams);
    depthI = shapeI(1,1);
    plot(shapeI(:,2), shapeI(:,1),'Color',cmap(1,:)); hold on;
    set(gca,'YDir','Reverse'); xlabel('Radius (m)'); ylabel('Depth (m)');
    plot(-shapeI(:,2), shapeI(:,1),'Color',cmap(1,:));
    plot([-shapeI(1,2) shapeI(1,2)],[depthI depthI],'Color',cmap(1,:));
    
    % Best 3 crater geometries
    % crater A
    shapeFA=shapes.shapeFA;
    depthFA=shapes.depthFA;
    plot(shapeFA(:,2), shapeFA(:,1),'Color',cmap(2,:));
    plot(-shapeFA(:,2), shapeFA(:,1),'Color',cmap(2,:));
    plot([-shapeFA(1,2) shapeFA(1,2)],[depthFA depthFA],'Color',cmap(2,:));
    labelA=sprintf('Optimal, misfit=%.4f', best5_misfits(1));
    % crater B
    shapeFB=shapes.shapeFB;
    depthFB=shapes.depthFB;
    plot(shapeFB(:,2), shapeFB(:,1),'Color',cmap(3,:));
    plot(-shapeFB(:,2), shapeFB(:,1),'Color',cmap(3,:));
    plot([-shapeFB(1,2) shapeFB(1,2)],[depthFB depthFB],'Color',cmap(3,:));
    labelB=sprintf('Second best, misfit=%.4f', best5_misfits(2));
    % crater C
    shapeFC=shapes.shapeFC;
    depthFC=shapes.depthFC;
    plot(shapeFC(:,2), shapeFC(:,1),'Color',cmap(4,:));
    plot(-shapeFC(:,2), shapeFC(:,1),'Color',cmap(4,:));
    plot([-shapeFC(1,2) shapeFC(1,2)],[depthFC depthFC],'Color',cmap(4,:));
    labelC=sprintf('Third best, misfit=%.4f', best5_misfits(3));
    
    legend('Initial','','',labelA,'','',labelB,'','',labelC);
    
    % Spectra
    figure(2);
    plot(dataF, abs(dataAmp)./max(abs(dataAmp)),'k','LineWidth',2);
    hold on; xlim([0 3]);
    xlabel('Frequency (Hz)');
    ylabel('Normalized Amplitude Spectra');
    plot(f, abs(spec.finA),'Color',cmap(2,:)); % crater A
    plot(f, abs(spec.finB),'Color',cmap(3,:)); % crater B
    plot(f, abs(spec.finC),'Color',cmap(4,:)); % crater B
    legend('Data','Optimal Model', 'Second best', 'Third best');
    
    % Misfit history
    figure(3);
    plot(misfit,'LineWidth',2);
    ylabel('Misfit');
    xlabel('Iteration Number');
    title('Misfit Evolution');
    grid on;
    
    % Parameter evolution
    figure(4);
    k = length(params);
    for i = 1:k
        subplot(k,1,i);
        plot(x(:,i) .* params(i),'LineWidth',2);
        ylabel(sprintf('Param %d',i));
        hold on;
        yline(lowerBnds(i),'r--');
        yline(upperBnds(i),'r--');
        grid on;
    end
    xlabel('Iteration Number');
end