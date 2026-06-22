% Plots SA data. Uses data saved from sim_anneal_data_inversion script
% Cole Campbell, 21st January 2026

clc, clearvars, close all;
cmap = get(gca,'ColorOrder');
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',18);

% load data
load("simulated_annealing_output_data\SAInvOut_Etna2018phase1_runlength600min_date2026-01-23_07_20_40_R0100m_D200m_T100C_freqLim3_nx5.mat");

%% process data
dataF = data_freq(:,1); % unpack data_freq
dataAmp = data_freq(:,2);


%% Print best basin datafprintf('Best %d distinct basins:\n', nBasins);
for i = 1:length(best5_misfits)
    fprintf('Basin %d: misfit = %.6f (local min at iteration %d), ', ...
        i, best5_misfits(i), best5_indices(i));
    fprintf('Parameters: [%s]\n', num2str(best5_params(i, :).*geomParams, '%.4f '));
end



%% plot figures - copied from sim_anneal_data_inversion
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
plot(dataF, abs(dataAmp)./max(abs(dataAmp)),'Color',cmap(1,:),'LineWidth',2);
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