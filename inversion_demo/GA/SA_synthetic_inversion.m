%% SA SYNTHETIC INVERSION %%
%
% Script file that performs SA inversion to invert harmonic infrasound
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
addpath ('C:\Users\dsp11\OneDrive\Summer Research\Github Files\inversion_demo')
addpath ('C:\Users\dsp11\OneDrive\Summer Research\Github Files\source\resonance')
addpath ('C:\Users\dsp11\OneDrive\Summer Research\Github Files\source\SBPoperators')
addpath ('C:\Users\dsp11\OneDrive\Summer Research\Github Files\source\inv')

%% Logicals %%
save_output = 1; % logical that determines if outputs are saved or not
plot_output = 1; % logical that determines if outputs are plotted
plot_input = 1; % logical that determines if synthetic inputs are plotted
rand_inputs = 0; % logical that determines if a random input is used

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

%% Define crater geometry (true parameters)
geomFlag = 1; % set for inversion
geomLowerBnds(1) = 50; % min depth
geomUpperBnds(1) = 200; % max depth
geomLowerBnds(2:8) = 30;   % minimum radius
geomUpperBnds(2:8) = 120;  % maximum radius
nGeom = 8;
if rand_inputs
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

%% Example: Call simulannealbnd for your problem

% Preallocate history
misfitHistory = [];

% Wrap the misfit calculation in a function that records history
fitnessFunc = @(x) trackMisfit(x);

function fval = trackMisfit(x)
    % Compute misfit for parameters x
    % Extract geom and source parameters from x
    if geomFlag
        geomParams_current = x(1:nParamsGeom);
    else
        geomParams_current = geomParams;
    end
    if srcFlag
        srcParams_current = x(end);
    else
        srcParams_current = srcParams;
    end

    try
        % Forward model evaluation (same as in GA)
        shape = geomFunction(geomParams_current);
        depth = shape(1,1);
        if depth <= 0 || any(shape(:,2) <= 0)
            fval = 1e10;
        else
            [S, ~, ~, ~] = sourceFunction(1, srcParams_current, srcStyle, discrParams);

            style = 'baffled piston';
            order = 4;
            res = resonance1d(shape, depth, [0 Nyquist], Nf, style, order, M);

            len_min = min(length(res.P), length(S));
            sim_P = res.P(1:len_min) .* S(1:len_min);

            sim_pFreq = [sim_P; conj(sim_P(end-1:-1:2))];
            sim_pFreq(length(sim_P)+1) = real(sim_pFreq(length(sim_P)+1));
            sim_pTime = ifft(sim_pFreq,'symmetric') / dt;

            sim_pTimeFilt = bandpass_butterworth(sim_pTime, filterband, Fs, filterorder);

            L = length(sim_pTimeFilt);
            NFFT = 2^nextpow2(L);
            sim_pFreqFilt = fft(sim_pTimeFilt,NFFT)/L;
            sim_fFilt = Fs/2 * linspace(0,1,NFFT/2+1);

            simNorm = abs(sim_pFreqFilt(1:NFFT/2+1));
            simNorm = simNorm / max(simNorm);

            % Interpolate to data frequencies
            simInterp = interp1(sim_fFilt, simNorm, dataF, 'linear', 0);

            idx = dataF <= freqLim;
            fval = sum((dataNorm(idx) - simInterp(idx)).^2);
        end
    catch
        fval = 1e10;
    end

    % Record misfit history
    misfitHistory(end+1,1) = fval;
end

% Define bounds and initial guess
x0 = geomParams;
lb = geomLowerBnds;
ub = geomUpperBnds;

% Create options for simulannealbnd
options = optimoptions('simulannealbnd', ...
    'Display', 'off', ...               % suppress output
    'MaxIterations', nGen, ...          % number of generations
    'FunctionTolerance', 1e-12);        % convergence tolerance

% Run Simulated Annealing
[xBest, fBest] = simulannealbnd(fitnessFunc, x0, lb, ub, options);

% Return outputs compatible with GA outputs
population = [];              % SA doesn't maintain a full population
misfit = misfitHistory;       % misfit per iteration
simSpec = generateSimSpec(xBest, srcParams, srcStyle, filterProps, data_freq); % optional
f = dataF;                     % frequency vector from your data
bestIndiv = xBest;             % best solution found

