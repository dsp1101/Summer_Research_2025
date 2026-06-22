%% GA SYNTHETIC INVERSION %%
%
% Script file that performs GA inversion to invert harmonic infrasound
% observations for crater geometry using synthetic data generated using
% CRes
%
% Written by Daniel Spencer (Adapted from Leighton Watson)
% December, 2025
% dsp65@uclive.ac.nz // daniel.spencer2007@gmail.com


%% Function setup %%
clear all; clc;
cmap = get(gca,'ColorOrder');
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',18);

% Find resonance, SBPoperators and inv folders
path(pathdef)
addpath inversion_demo\GA\
addpath ('C:\Users\dsp11\OneDrive\Summer Research 25-26\Github Files\inversion_demo')
addpath ('C:\Users\dsp11\OneDrive\Summer Research 25-26\Github Files\source\resonance')
addpath ('C:\Users\dsp11\OneDrive\Summer Research 25-26\Github Files\source\SBPoperators')
addpath ('C:\Users\dsp11\OneDrive\Summer Research 25-26\Github Files\source\inv')

%% Logicals %%
save_output = 1; % logical that determines if outputs are saved or not
plot_output = 1; % logical that determines if outputs are plotted
plot_input = 1; % logical that determines if synthetic inputs are plotted
rand_synth_geom = 0; % logical that determines if a random input is used

%% Discretization parameters

% set the parameters for the resonance1d calculations 
T = 25; % total time (s)
N = 250; % number of grid points (formulas assume even N)
dt = T/N; % time step (s)
Nyquist = 1/(2*dt); % Nyquist frequency (Hz)
Nf = N/2+1; % number of frequency samples
discrParams = [T N Nf Nyquist dt]; % save parameters into array

%% Temperature parameters

craterTemp = 100; % crater temperature
atmoTemp = 0; % atmospheric temperature
temp = [craterTemp, atmoTemp]; 
M = problemParametersInv(craterTemp, atmoTemp);

%% GA Parameters
useTimeLimit = true;
maxTime = 60;
popSize = 10;          % population size
nGen = 5;            % number of generations
mutationRate = 0.1;    % probability of mutation
crossoverRate = 0.8;   % probability of crossover
eliteCount = 2;        % number of elite individuals to preserve
tournamentSize = 3;    % tournament selection size

%% Define crater geometry (true parameters)
geomFlag = 1; % set for inversion
geomLowerBnds(1) = 50; % min depth
geomUpperBnds(1) = 200; % max depth
geomLowerBnds(2:8) = 30;   % minimum radius
geomUpperBnds(2:8) = 120;  % maximum radius
nGeom = 8;
if rand_synth_geom
    geomParams = geomLowerBnds + rand(1, nGeom) .* (geomUpperBnds - geomLowerBnds); % Generate random inputs within set bounds
else
    geomParams = [150 50 80 80 80 80 80 90]; % manualy set inputs [depth R0 R1 ... R7]
    % note: R0 is radius at max depth, R7 is radius at zero depth
end
geomR0 = geomParams(2);
geomDepth = geomParams(1);
shape = geomFunction(geomParams);
depth = shape(1,1);

%% Define Source
srcFlag = 0; % dont invert source
srcParams = 0.2;
srcStyle = 'Gauss';
[S, ~, ~, ~] = sourceFunction(1, srcParams, srcStyle, discrParams);

%% Compute resonance transfer function
style = 'baffled piston';
order = 4;
freq = [0 Nyquist];
res = resonance1d(shape, depth, freq, Nf, style, order, M);

%% Convolve source and resonance
sim_P = res.P(1:N/2+1) .* S(1:N/2+1);
sim_P = sim_P(:);

%% Convert to time domain
sim_pFreq = [sim_P; conj(sim_P(end-1:-1:2))];
sim_pFreq(N/2+1) = real(sim_pFreq(N/2+1));
sim_pTime = ifft(sim_pFreq, 'symmetric') / dt;

%% Apply filtering (bandpass)
filterband = [0.01 3]; % Hz
filterorder = 4;
Fs = 1/dt; % correct sampling rate
filterProps = [filterband filterorder Fs];
freqLim = 3;
sim_pTimeFilt = bandpass_butterworth(sim_pTime, filterband, Fs, filterorder);

%% Convert back to frequency domain
L = length(sim_pTimeFilt);
NFFT = 2^nextpow2(L);
sim_pFreqFilt = fft(sim_pTimeFilt, NFFT) / L;
dataF = Fs/2 * linspace(0, 1, NFFT/2+1)';
dataAmp = abs(sim_pFreqFilt(1:NFFT/2+1));
dataAmp = dataAmp / max(dataAmp); % normalize

%% Combine into single variable
data_freq = [dataF, dataAmp];

%% run GA inversion %%
tic
[population, misfit, simSpec, f, bestIndiv] = ...
    ga_spec(popSize, nGen, mutationRate, crossoverRate, eliteCount, ...
            tournamentSize, geomParams, geomFlag, srcParams, srcFlag, ...
            srcStyle, geomLowerBnds, geomUpperBnds, discrParams, temp, ...
            filterProps, data_freq, freqLim, ...
            useTimeLimit, maxTime);
toc

%% Compute best geometry %%
bestParams = geomLowerBnds + bestIndiv(end,:) .* (geomUpperBnds - geomLowerBnds);
shapeF = geomFunction(bestParams);
depthF = shapeF(1,1);

% Synthetic crater geometry (true parameters used to generate data)
syntheticShape = geomFunction(geomParams);
syntheticDepth = syntheticShape(1,1);

%% Plot outputs %%
if plot_output
    % 1. Crater geometry (synthetic vs GA best-fit)
    figure(1); hold on;
    % Synthetic crater in light blue
    fill([syntheticShape(:,2); -flipud(syntheticShape(:,2))], ...
         [syntheticShape(:,1); flipud(syntheticShape(:,1))], ...
         [0.6 0.8 1], 'FaceAlpha',0.5,'EdgeColor','b','LineWidth',1.5);
    % GA best-fit crater in red outline
    plot(shapeF(:,2), shapeF(:,1),'r','LineWidth',2);
    plot(-shapeF(:,2), shapeF(:,1),'r','LineWidth',2);
    % Top edge
    plot([-shapeF(1,2) shapeF(1,2)],[depthF depthF],'r','LineWidth',2);

    xlabel('Radius (m)'); ylabel('Depth (m)');
    title('Crater Geometry: Synthetic (blue) vs GA Best Fit (red)');
    axis equal; grid on;

    % 2. Spectrum comparison (synthetic vs GA best fit)
    figure(2); hold on;
    plot(dataF, dataAmp,'b','LineWidth',1.5);   % synthetic data in blue
    plot(f, abs(simSpec),'r','LineWidth',2);    % GA best fit in red
    xlabel('Frequency (Hz)'); ylabel('Normalized Amplitude');
    title('Synthetic Spectrum vs GA Best Fit');
    legend('Synthetic Data','GA Best Fit'); grid on;

    % 3. Misfit evolution
    figure(3); clf;
    plot(misfit,'r' ,'LineWidth',2); xlabel('Generation'); ylabel('Best Misfit'); 
    title('GA Misfit Evolution'); grid on;
end
%% Save outputs %%
dataStr = 'synthetic_test1';
if save_output
    pathname = 'C:\Users\dsp11\OneDrive\Summer Research\Saved Results GA'; % Destination folder of results
    filename = ['GA_InvOut_',dataStr,'_R0',num2str(geomR0),'_D',num2str(geomDepth)];
    save(fullfile(pathname,filename),'population','misfit','simSpec','f','bestIndiv','bestParams');
end