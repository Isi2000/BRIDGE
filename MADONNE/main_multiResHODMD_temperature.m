%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%         Multi-dimensional HODMD Analysis for Temperature Data       %%%
%%%      Extension of HODMD for the analysis of temperature fields      %%%
%%%                   Adapted for binned_temperature_cs.mat             %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Adapted from the original HODMD implementation by Le Clainche & Vega   %
%  for temperature field analysis instead of flow field analysis          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
%%%%%%%%%%%%%%%%%%%%%% SAVE RESULTS IN FOLDER DMD_solution %%%%%%%%%%%%%%%
mkdir('DMD_solution_temperature')
system('rm -r DMD_solution_temperature');
mkdir('DMD_solution_temperature')
filename = sprintf('./DMD_solution_temperature/DMD_history.txt' );
diary(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% USER:
% Load temperature data in tensor form: binned_temperature_cs
% Dimensions: 41x76x161x500
% i: X spatial dimension (41 points)
% j: Y spatial dimension (76 points) 
% k: Z spatial dimension (161 points)
% t: Time dimension (500 time steps)

load /home/isacco/DATASETS/BINNED_TEMPERATURES/binned_temperature_cs.mat
Tensor = binned_temperature_cs;

%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the snapshot number (using all available time steps)
SNAP = size(Tensor, 4);  % 500 time steps
varepsilon1 = 5e-2;      % First tolerance (SVD)
varepsilon2 = 5e-2;      % Second tolerance (DMD-d modes)
d = 100;                 % Parameter of DMD-d (higher order Koopman assumption)
deltaT = 1.0;            % Time between snapshots (adjust as needed)

% Set the position of the temporal variable. For temperature data,
% time is in the 4th dimension
TimePos = 4;

fprintf('Original tensor dimensions: %s\n', mat2str(size(Tensor)));
fprintf('Number of snapshots: %d\n', SNAP);
fprintf('Time step: %f\n', deltaT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Tensor dimension - number of snapshots: 
Tensor0 = Tensor;
clear Tensor

% For 4D tensor (3 spatial + 1 temporal)
if TimePos == 4
    Tensor(:,:,:,:) = Tensor0(:,:,:,1:SNAP);
elseif TimePos == 5
    Tensor(:,:,:,:,:) = Tensor0(:,:,:,:,1:SNAP);
end

Time = [1:SNAP] * deltaT;

fprintf('Processing tensor with dimensions: %s\n', mat2str(size(Tensor)));

%% ALGORITHM:
%% ITERATIVE
nn0 = size(Tensor);
nn(1) = nn0(1);
nn(2:length(nn0)) = 0;

fprintf('Starting HOSVD decomposition...\n');
%% PERFORM HOSVD DECOMPOSITION TO CALCULATE THE REDUCED TEMPORAL MATRIX hatT
[hatT, U, S, sv, nn1] = hosvd_function(Tensor, varepsilon1, nn, nn0, TimePos);
% hatT: reduced temporal matrix
% U: temporal SVD modes
% S: tensor core
% sv: singular values
% nn: number of singular values

fprintf('HOSVD completed. Reduced temporal matrix size: %s\n', mat2str(size(hatT)));

%% PERFORM DMD-d TO THE REDUCED TEMPORAL MATRIX hatT
fprintf('Starting DMD analysis (d = %d)...\n', d);
if d > 1
    [GrowthRate, Frequency, Amplitude, hatMode] = DMDd(d, hatT, Time, varepsilon1, varepsilon2);
else
    [GrowthRate, Frequency, Amplitude, hatMode] = DMD1(hatT, Time, varepsilon1, varepsilon2);
end

fprintf('DMD analysis completed. Found %d modes.\n', length(Frequency));

%% RECONSTRUCT THE ORIGINAL TENSOR USING THE DMD EXPANSION
fprintf('Reconstructing temperature field...\n');
[TensorReconst] = DMDreconst(GrowthRate, Frequency, hatMode, Time, U, S, sv, nn1, TimePos);

% For temperature data, take real part
TensorReconst = real(TensorReconst);

%% CALCULATE ERROR METRICS
fprintf('Calculating reconstruction error...\n');
RRMSE = norm(Tensor(:) - TensorReconst(:), 2) / norm(Tensor(:), 2);
fprintf('Relative mean square error: %.6e\n', RRMSE);

% Display results
fprintf('\n=== TEMPERATURE FIELD DMD ANALYSIS RESULTS ===\n');
fprintf('Growth rate, Frequency, Amplitude:\n');
GrowthrateFrequencyAmplitude = [GrowthRate', Frequency', Amplitude'];
for i = 1:min(10, length(Frequency))  % Show first 10 modes
    fprintf('Mode %2d: δ=%.4f, ω=%.4f, a=%.4e\n', i, GrowthRate(i), Frequency(i), Amplitude(i));
end
if length(Frequency) > 10
    fprintf('... and %d more modes\n', length(Frequency) - 10);
end

diary off

%% SAVE RESULTS
fprintf('Saving results to DMD_solution_temperature/...\n');
% Save the reconstruction and analysis results
save ./DMD_solution_temperature/TensorReconst_temperature.mat TensorReconst -v7.3
save ./DMD_solution_temperature/GrowthrateFrequencyAmplitude_temperature.mat GrowthrateFrequencyAmplitude

fprintf('Calculating temperature DMD modes...\n');
%% Calculate DMD modes
[N, ~] = size(hatT);
[DMDmode] = calculateDMDmode(N, hatMode, Amplitude, U, S, nn1, TimePos);
% Save DMD modes
save ./DMD_solution_temperature/DMDmode_temperature.mat DMDmode -v7.3

%% VISUALIZATION
fprintf('Creating visualizations...\n');

% Frequency vs. absolute value of Growth rate
h1 = figure('Name', 'Temperature Field: Frequency vs Growth Rate');
semilogy(Frequency, abs(GrowthRate), 'o', 'linewidth', 2, 'color', 'r');
grid on;
title('Temperature Field: DMD Spectrum');
xlabel('Frequency \omega_m');
ylabel('|Growth Rate \delta_m|');
name1 = sprintf('./DMD_solution_temperature/FrequencyGrowthrate_temperature_d%03i', d);
saveas(h1, name1, 'fig');
saveas(h1, [name1 '.png']);

% Frequency vs. amplitudes
h2 = figure('Name', 'Temperature Field: Frequency vs Amplitude');
semilogy(Frequency, Amplitude, 'o', 'linewidth', 2, 'color', 'b');
grid on;
title('Temperature Field: DMD Mode Amplitudes');
xlabel('Frequency \omega_m');
ylabel('Amplitude a_m');
name2 = sprintf('./DMD_solution_temperature/FrequencyAmplitude_temperature_d%03i', d);
saveas(h2, name2, 'fig');
saveas(h2, [name2 '.png']);

% Time evolution of first few modes
h3 = figure('Name', 'Temperature Field: Temporal Evolution');
t_plot = Time(1:min(100, length(Time)));  % Plot first 100 time steps
hold on;
colors = {'r', 'g', 'b', 'm', 'c'};
for i = 1:min(5, length(Frequency))
    mode_evolution = real(Amplitude(i) * exp((GrowthRate(i) + 1i*Frequency(i)) * t_plot));
    plot(t_plot, mode_evolution, colors{i}, 'linewidth', 1.5, ...
         'DisplayName', sprintf('Mode %d (ω=%.3f)', i, Frequency(i)));
end
grid on;
title('Temperature Field: Temporal Evolution of DMD Modes');
xlabel('Time');
ylabel('Mode Amplitude');
legend('Location', 'best');
name3 = sprintf('./DMD_solution_temperature/TemporalEvolution_temperature_d%03i', d);
saveas(h3, name3, 'fig');
saveas(h3, [name3 '.png']);

%% SUMMARY STATISTICS
fprintf('\n=== SUMMARY STATISTICS ===\n');
fprintf('Dataset: binned_temperature_cs.mat\n');
fprintf('Original dimensions: %s\n', mat2str(size(binned_temperature_cs)));
fprintf('Processed dimensions: %s\n', mat2str(size(Tensor)));
fprintf('Number of DMD modes found: %d\n', length(Frequency));
fprintf('Reconstruction error (RRMSE): %.6e\n', RRMSE);
fprintf('Dominant frequency: %.4f\n', Frequency(1));
fprintf('Dominant growth rate: %.4f\n', GrowthRate(1));

% Find stable, growing, and decaying modes
stable_modes = sum(abs(GrowthRate) < 0.01);
growing_modes = sum(GrowthRate > 0.01);
decaying_modes = sum(GrowthRate < -0.01);

fprintf('\nMode Classification:\n');
fprintf('  Stable modes (|δ| < 0.01): %d\n', stable_modes);
fprintf('  Growing modes (δ > 0.01): %d\n', growing_modes);
fprintf('  Decaying modes (δ < -0.01): %d\n', decaying_modes);

fprintf('\nAnalysis completed successfully!\n');
fprintf('Check the DMD_solution_temperature/ folder for all results.\n');