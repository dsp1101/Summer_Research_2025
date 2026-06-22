%% Graph Maker %%
%% ================================
% INVERSION RESULTS PLOTTING (GA/MCMC/SA/Mixed)
% ================================

clear; clc; close all;
cd("C:\Users\dsp11\OneDrive\Summer Research 25-26\Saved Results GA\Data Inversion Results")

%% ---- USER CONFIGURATION ----
% Choose which dataset to load:
dataset_type = 'combined';  % Options: 'GA', 'MCMC', 'combined'

% For combined dataset, specify file name
if strcmp(dataset_type, 'combined')
    filename = 'all_results_inversion.csv';  % Your combined MCMC/GA/SA/Mixed file
else
    filename = 'Data Inversion Results MATH395 25-26.csv';  % Your GA-only file
end

%% ---- LOAD DATA ----
data = readtable(filename);

% Check if this is the combined dataset (has 'Model' column)
if ismember('Model', data.Properties.VariableNames)
    % Combined dataset (MCMC/GA/SA/Mixed)
    data.model = data.Model;
    data.minutes = data.Minutes;
    data.dataset = data.Dataset;
    data.best_misfit = data.BestMisfit;
    data.best_geometry = data.BestGeometry;
    
    % For combined data, we don't have GA-specific parameters
    has_ga_params = false;
else
    % GA-only dataset
    data(data.Minutes == 600, :) = [];
    data.model = repmat({'GA'}, height(data), 1);  % Add model column
    data.minutes = data.Minutes;
    data.dataset = data.Dataset;
    data.best_misfit = data.BestMisfit;
    data.generations_completed = data.GenerationsCompleted;
    data.best_geometry = data.BestGeometry;
    data.population_size = data.PopSize;
    data.mutation_rate = data.MutationRate;
    data.crossover = data.CrossoverRate;
    data.elite_count = data.EliteCount;
    data.tournament_size = data.TournamentSize;
    
    has_ga_params = true;
end

disp('Available columns:');
disp(data.Properties.VariableNames);
fprintf('\nLoaded %d entries\n', height(data));

% Get unique models
models = unique(data.model);
fprintf('Models in dataset: ');
fprintf('%s ', models{:});
fprintf('\n\n');

%% ================================
% 1) MISFIT vs RUNTIME (COLOURED BY MODEL)
%% ================================
figure('Position', [100, 100, 900, 600]);
hold on;

colors_model = lines(length(models));

for i = 1:length(models)
    idx = strcmp(data.model, models{i});
    scatter(data.minutes(idx), ...
            data.best_misfit(idx), ...
            100, colors_model(i,:), 'filled', ...
            'DisplayName', models{i});
end

xlabel('Runtime (minutes)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Best Misfit', 'FontWeight', 'bold', 'FontSize', 12);
title('Best Misfit vs Runtime (by Model)', 'FontWeight', 'bold', 'FontSize', 14);
legend('Location', 'best');
grid on;
box on;

%% ================================
% 2) MISFIT vs RUNTIME (COLOURED BY DATASET)
%% ================================
figure('Position', [100, 100, 900, 600]);
hold on;

datasets = unique(data.dataset);
colors_dataset = lines(length(datasets));

for i = 1:length(datasets)
    idx = strcmp(data.dataset, datasets{i});
    scatter(data.minutes(idx), ...
            data.best_misfit(idx), ...
            100, colors_dataset(i,:), 'filled', ...
            'DisplayName', datasets{i});
end

xlabel('Runtime (minutes)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Best Misfit', 'FontWeight', 'bold', 'FontSize', 12);
title('Best Misfit vs Runtime (by Dataset)', 'FontWeight', 'bold', 'FontSize', 14);
legend('Location', 'best', 'Interpreter', 'none');
grid on;
box on;

%% ================================
% 3) MODEL COMPARISON BOXPLOTS
%% ================================
figure('Position', [100, 100, 900, 600]);

subplot(1,2,1);
boxplot(data.best_misfit, data.model);
ylabel('Best Misfit', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Model', 'FontWeight', 'bold', 'FontSize', 12);
title('Misfit Distribution by Model', 'FontWeight', 'bold', 'FontSize', 14);
grid on;
box on;

subplot(1,2,2);
boxplot(data.best_misfit, data.dataset);
ylabel('Best Misfit', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Dataset', 'FontWeight', 'bold', 'FontSize', 12);
title('Misfit Distribution by Dataset', 'FontWeight', 'bold', 'FontSize', 14);
grid on;
box on;
xtickangle(45);

%% ================================
% 4) VIOLIN/STRIP PLOT: MISFIT BY MODEL AND DATASET
%% ================================
figure('Position', [100, 100, 1200, 600]);

% Create combinations
subplot(1,2,1);
hold on;
model_colors = lines(length(models));

for i = 1:length(models)
    for j = 1:length(datasets)
        idx = strcmp(data.model, models{i}) & strcmp(data.dataset, datasets{j});
        if sum(idx) > 0
            x_pos = (i-1)*length(datasets) + j;
            scatter(repmat(x_pos, sum(idx), 1), data.best_misfit(idx), ...
                    80, model_colors(i,:), 'filled', 'MarkerFaceAlpha', 0.6);
        end
    end
end

xticks(1:length(models)*length(datasets));
xticklabels_custom = {};
for i = 1:length(models)
    for j = 1:length(datasets)
        xticklabels_custom{end+1} = sprintf('%s\n%s', models{i}, datasets{j});
    end
end
xticklabels(xticklabels_custom);
xtickangle(45);
ylabel('Best Misfit', 'FontWeight', 'bold', 'FontSize', 12);
title('Misfit by Model-Dataset Combination', 'FontWeight', 'bold', 'FontSize', 14);
grid on;
box on;

% Runtime comparison
subplot(1,2,2);
hold on;
for i = 1:length(models)
    idx = strcmp(data.model, models{i});
    scatter(data.minutes(idx), data.best_misfit(idx), ...
            100, model_colors(i,:), 'filled', ...
            'DisplayName', models{i});
end
xlabel('Runtime (minutes)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Best Misfit', 'FontWeight', 'bold', 'FontSize', 12);
title('Model Performance vs Runtime', 'FontWeight', 'bold', 'FontSize', 14);
legend('Location', 'best');
grid on;
box on;

%% ================================
% 5) GA-SPECIFIC PLOTS (only if GA parameters available)
%% ================================
if has_ga_params
    % Generations vs Runtime
    figure('Position', [100, 100, 800, 600]);
    scatter(data.minutes, data.generations_completed, 80, 'filled');
    xlabel('Runtime (minutes)', 'FontWeight', 'bold', 'FontSize', 12);
    ylabel('Generations Completed', 'FontWeight', 'bold', 'FontSize', 12);
    title('Generations vs Runtime', 'FontWeight', 'bold', 'FontSize', 14);
    grid on;
    box on;
    
    % Parameter sweeps
    figure('Position', [100, 100, 1400, 900]);
    
    subplot(2,3,1);
    scatter(data.population_size, data.best_misfit, 80, 'filled');
    xlabel('Population Size', 'FontWeight', 'bold');
    ylabel('Best Misfit', 'FontWeight', 'bold');
    title('Population Size vs Misfit');
    grid on; box on;
    
    subplot(2,3,2);
    scatter(data.mutation_rate, data.best_misfit, 80, 'filled');
    xlabel('Mutation Rate', 'FontWeight', 'bold');
    ylabel('Best Misfit', 'FontWeight', 'bold');
    title('Mutation Rate vs Misfit');
    grid on; box on;
    
    subplot(2,3,3);
    scatter(data.crossover, data.best_misfit, 80, 'filled');
    xlabel('Crossover Rate', 'FontWeight', 'bold');
    ylabel('Best Misfit', 'FontWeight', 'bold');
    title('Crossover Rate vs Misfit');
    grid on; box on;
    
    subplot(2,3,4);
    scatter(data.elite_count, data.best_misfit, 80, 'filled');
    xlabel('Elite Count', 'FontWeight', 'bold');
    ylabel('Best Misfit', 'FontWeight', 'bold');
    title('Elite Count vs Misfit');
    grid on; box on;
    
    subplot(2,3,5);
    scatter(data.tournament_size, data.best_misfit, 80, 'filled');
    xlabel('Tournament Size', 'FontWeight', 'bold');
    ylabel('Best Misfit', 'FontWeight', 'bold');
    title('Tournament Size vs Misfit');
    grid on; box on;
    
    subplot(2,3,6);
    scatter(data.generations_completed, data.best_misfit, 80, 'filled');
    xlabel('Generations', 'FontWeight', 'bold');
    ylabel('Best Misfit', 'FontWeight', 'bold');
    title('Generations vs Misfit');
    grid on; box on;
    
    % Correlation heatmap
    vars = data{:,{'minutes','best_misfit','generations_completed', ...
                   'population_size','mutation_rate','crossover', ...
                   'elite_count','tournament_size'}};
    
    labels = {'Minutes','Misfit','Generations','Pop','Mut','Cross','Elite','Tourn'};
    
    C = corrcoef(vars,'Rows','complete');
    
    figure('Position', [100, 100, 700, 600]);
    imagesc(C);
    colorbar;
    colormap('jet');
    caxis([-1 1]);
    axis equal tight;
    xticks(1:length(labels));
    yticks(1:length(labels));
    xticklabels(labels);
    yticklabels(labels);
    xtickangle(45);
    title('Correlation Matrix (GA Parameters)', 'FontWeight', 'bold', 'FontSize', 14);
    
    % Add correlation values as text
    for i = 1:length(labels)
        for j = 1:length(labels)
            text(j, i, sprintf('%.2f', C(i,j)), ...
                'HorizontalAlignment', 'center', ...
                'Color', 'w', 'FontWeight', 'bold');
        end
    end
end

%% ================================
% PLOT SPECIFIC GEOMETRY AND SPECTRA
%% ================================

% USER INPUT: Specify which entry to plot
entry_to_plot = 22;  % Change this to plot different entries

fprintf('\n=== Plotting Entry %d ===\n', entry_to_plot);
fprintf('Model: %s\n', data.model{entry_to_plot});
fprintf('Dataset: %s\n', data.dataset{entry_to_plot});
fprintf('Misfit: %.4f\n', data.best_misfit(entry_to_plot));
fprintf('Runtime: %.2f min\n', data.minutes(entry_to_plot));

% Parse geometry string (handles both "[depth, r1, ...]" and "depth r1 ..." formats)
geom_str = data.best_geometry{entry_to_plot};

% Remove brackets if present
geom_str = strrep(geom_str, '[', '');
geom_str = strrep(geom_str, ']', '');

% Remove any extra whitespace
geom_str = strtrim(geom_str);

% Try to parse as numbers
best_params = str2num(geom_str);

% If that fails, try splitting by spaces or commas
if isempty(best_params)
    % Try comma-separated
    parts = strsplit(geom_str, ',');
    if length(parts) == 1
        % Try space-separated
        parts = strsplit(geom_str);
    end
    
    % Convert each part to number
    best_params = [];
    for p = 1:length(parts)
        num = str2double(strtrim(parts{p}));
        if ~isnan(num)
            best_params(end+1) = num;
        end
    end
end

if isempty(best_params) || length(best_params) < 6
    warning('Could not parse geometry for entry %d (got %d parameters)', ...
            entry_to_plot, length(best_params));
    fprintf('Raw string: "%s"\n', data.best_geometry{entry_to_plot});
else
    fprintf('Geometry: ');
    fprintf('%.2f ', best_params);
    fprintf('\n');
    
    % Load the corresponding observational data
    dataset_name = data.dataset{entry_to_plot};
    
    % Add path to your data files
    data_path = 'C:\Users\dsp11\OneDrive\Summer Research 25-26\Github Files\GRL2020\data';
    addpath(data_path);
    
    % Load observational data
    try
        datafile = strcat(dataset_name, '.mat');
        load(fullfile(data_path, datafile));
        
        % Run forward model to get simulated spectrum
        % (You'll need to adjust these parameters to match your setup)
        craterTemp = 100;
        atmoTemp = 0;
        temp = [craterTemp, atmoTemp];
        
        T = 25;
        N = 250;
        dt = T/N;
        Nyquist = 1/(2*dt);
        Nf = N/2+1;
        freq = [0 Nyquist];
        discrParams = [T N Nf Nyquist dt];
        
        filterband = [0.25 4.8];
        filterorder = 4;
        Fs = 10;
        
        order = 4;
        style = 'baffled piston';
        M = problemParametersInv(craterTemp, atmoTemp);
        
        % Generate geometry
        shape = geomFunction(best_params);
        depth = shape(1,1);
        
        % Source function
        srcParams = 0.3;
        srcStyle = 'Brune';
        [S, ~, ~, ~] = sourceFunction(1, srcParams, srcStyle, discrParams);
        
        % Compute transfer function
        res = resonance1d(shape, depth, freq, Nf, style, order, M);
        
        % Convolve and normalize
        sim.f = res.f;
        sim.P = (res.P(1:N/2+1).*S(1:N/2+1))./max(res.P(1:N/2+1).*S(1:N/2+1));
        
        % Time domain processing
        sim.pFreq = [sim.P conj(sim.P(end-1:-1:2))];
        sim.pFreq(N/2+1) = real(sim.pFreq(N/2+1));
        sim.pTime = ifft(sim.pFreq, 'symmetric')/dt;
        
        % Filter
        sim.pTimeFilt = bandpass_butterworth(sim.pTime, filterband, Fs, filterorder);
        
        % Back to frequency
        L = length(sim.pTimeFilt);
        NFFT = 2^nextpow2(L);
        sim.pFreqFilt = fft(sim.pTimeFilt, NFFT) / L;
        f_sim = Fs/2 * linspace(0, 1, NFFT/2+1);
        simSpec = abs(sim.pFreqFilt(1:NFFT/2+1));
        simSpec = simSpec / max(simSpec);
        
        % Plot Geometry
        figure('Position', [100, 100, 800, 600]);
        hold on;
        
        % Plot crater shape
        plot(shape(:,2), shape(:,1), 'Color', [0.1, 0.7, 0.1], 'LineWidth', 3);
        plot(-shape(:,2), shape(:,1), 'Color', [0.1, 0.7, 0.1], 'LineWidth', 3);
        plot([-shape(1,2) shape(1,2)], [depth depth], 'Color', [0.1, 0.7, 0.1], 'LineWidth', 3);
        
        set(gca, 'YDir', 'Reverse');
        xlabel('Radius (m)', 'FontWeight', 'bold', 'FontSize', 12);
        ylabel('Depth (m)', 'FontWeight', 'bold', 'FontSize', 12);
        title(sprintf('Crater Geometry - %s (Misfit: %.4f)', dataset_name, data.best_misfit(entry_to_plot)), ...
              'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'none');
        axis equal;
        grid on;
        box on;
        
        % Plot Spectra
        figure('Position', [100, 100, 800, 600]);
        hold on;
        
        % Observed data in blue
        plot(dataF, abs(dataAmp)./max(abs(dataAmp)), ...
             'Color', [0, 0.4470, 0.7410], 'LineWidth', 3);
        
        % Simulated in green
        plot(f_sim, simSpec, 'Color', [0.1, 0.7, 0.1], 'LineWidth', 2.5);
        
        xlim([0 3]);
        xlabel('Frequency (Hz)', 'FontWeight', 'bold', 'FontSize', 12);
        ylabel('Normalized Amplitude Spectrum', 'FontWeight', 'bold', 'FontSize', 12);
        title(sprintf('Observed vs Simulated - %s (Misfit: %.4f)', dataset_name, data.best_misfit(entry_to_plot)), ...
              'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'none');
        legend({'Observed Data', 'Simulated (GA Best Fit)'}, 'Location', 'best');
        grid on;
        box on;
        
    catch ME
        warning('Could not load or process data for dataset %s: %s', dataset_name, ME.message);
    end
end

%% ================================
disp("Plotting complete.");