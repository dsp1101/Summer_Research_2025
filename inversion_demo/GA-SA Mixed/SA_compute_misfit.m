function misfit = SA_compute_misfit(x_norm, params, geomParams, ...
    geomFlag, srcParams, srcFlag, srcStyle, discrParams, temp, ...
    filterProps, data_freq, freqLim)
% computes misfit value, for use in sim_anneal_data_inversion
% this is the objectiveFnc for simulannealbnd, a function in the Global
% Optimisation toolbox

vars = x_norm.*params;

%% CONSTRAINT CHECK: Monotonic decreasing radii %%
if geomFlag && length(geomParams) > 1
    geomLength = length(geomParams);
    geomParams_current = vars(1:geomLength);
    radii = geomParams_current(2:end);  % Extract radii (skip depth)
    
    % Check if radii are monotonically decreasing
    if any(diff(radii) > 0.01)  % Allow tiny tolerance for numerical errors
        misfit = 1e10;  % Return very large misfit for invalid geometry
        return;
    end
    
    % Also check for non-positive values
    if any(radii <= 0) || geomParams_current(1) <= 0
        misfit = 1e10;
        return;
    end
end

%% problem parameters %%
geomLength = length(geomParams);
srcLength = length(srcParams);
craterTemp = temp(1); % crater temperature
atmoTemp = temp(2); % atmosphere temperature
N = discrParams(2); % number of grid points
Nf = discrParams(3); % number of frequency samples
Nyquist = discrParams(4); % Nyquist frequency
freq = [0 Nyquist]; % frequency vector (Hz)
dt = discrParams(5); % time interval (s)
order = 4; % order of numerical scheme (4, 6 or 8)
style = 'baffled piston'; % acoustic radiation model ('monopole' or 'baffled piston')
M = problemParametersInv(craterTemp, atmoTemp); % problem parameters required for resonance1d

filterband = filterProps(1:2);
filterorder = filterProps(3);
Fs = filterProps(4); % sampling frequency

%% compute geometry
if geomFlag
    shape = geomFunction(vars(1:geomLength));
else
    shape = geomFunction(geomParams);
end
depth = shape(1,1);

%% compute source
if srcFlag
    if geomFlag
        [S, ~, ~, ~] = sourceFunction(1, vars(geomLength+1:geomLength+srcLength), srcStyle, discrParams);
    else
        [S, ~, ~, ~] = sourceFunction(1, vars(1:srcLength), srcStyle, discrParams);
    end
else
    [S, ~, ~, ~] = sourceFunction(1, srcParams, srcStyle, discrParams);
end

%% simulate spectra
try
    res = resonance1d(shape, depth, freq, Nf, style, order, M); % compute transfer function
    sim.f = res.f; % frequency vector
    sim.P = (res.P(1:N/2+1).*S(1:N/2+1))./max(res.P(1:N/2+1).*S(1:N/2+1)); % convolve transfer function and source function and normalize amplitudes
    
    % convert to time domain
    sim.pFreq = [sim.P conj(sim.P(end-1:-1:2))];
    sim.pFreq(N/2+1) = real(sim.pFreq(N/2+1));
    sim.pTime = ifft(sim.pFreq, 'symmetric')/dt;
    
    % apply filtering
    sim.pTimeFilt = bandpass_butterworth(sim.pTime, filterband, Fs, filterorder);
    
    % convert back to frequency domain
    L = length(sim.pTimeFilt);
    NFFT = 2^nextpow2(L);
    sim.pFreqFilt = fft(sim.pTimeFilt, NFFT)/L;
    sim.pFreqAbsNorm = abs(sim.pFreqFilt(1:NFFT/2+1))./max(abs(sim.pFreqFilt(1:NFFT/2+1)));
    sim.fFilt = Fs/2*linspace(0,1,NFFT/2+1);
    
    % data
    dataSpectrum.f = data_freq(:,1);
    dataSpectrum.P = abs(data_freq(:,2))./max(abs(data_freq(:,2)));
    
    % interpolate simulation onto same frequency vector as data
    sim.pFreqIntp = pchip(sim.fFilt, sim.pFreqAbsNorm, dataSpectrum.f);
    sim.F = dataSpectrum.f;
    
    % band limit signal
    idx = find(sim.F > freqLim, 1, 'first');
    if isempty(idx)
        idx = length(sim.F);
    end
    
    dataAmpSpec = abs(dataSpectrum.P(1:idx));
    simAmpSpec = abs(sim.pFreqIntp(1:idx));
    
    % Validate before computing misfit
    if any(isnan(simAmpSpec)) || any(isinf(simAmpSpec))
        misfit = 1e10;
        return;
    end
    
    % compute misfit
    misfit = norm(dataAmpSpec - simAmpSpec);
    
    % Final validation
    if isnan(misfit) || isinf(misfit)
        misfit = 1e10;
    end
    
catch ME
    % If forward model fails, return large misfit
    misfit = 1e10;
end
end