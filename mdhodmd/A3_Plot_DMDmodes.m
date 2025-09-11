%% CLEAR
clc; clear all; close all;

%% CALIBRATION
d_value = 50; % Number of HODMD windows
tol_value = 1e-3;      % Tolerance to truncate the number of modes retained

%% PLOT

% Load mode tensor and frequencies
load(sprintf('./Cases/mdHODMDit_solution_d%d_tol%0.1e/DMDmode.mat',d_value,tol_value));
load(sprintf('./Cases/mdHODMDit_solution_d%d_tol%0.1e/GrowthrateFrequencyAmplitude.mat',d_value,tol_value));

TimePos = ndims(DMDmode);
data_shape = size(DMDmode);

% Define the frequencies desired
freq = [0, 0.8, 1.6];
tolerance = 0.1; % Deviation permitted

% Find indices for desired frequencies
indices = NaN(size(freq));
for i = 1:length(freq)
    [min_val, idx] = min(abs(GrowthRateFrequencyAmplitude(:,2) - freq(i)));
    if min_val < tolerance
        indices(i) = idx;
    end
end

% Filter indices
if length(indices)>length(indices(~isnan(indices)))
    error('Some frequencies don´t match within tolerance. Check your input data.');
end
indices = indices(~isnan(indices));

% Plot modes
field = 1;

if TimePos == 4
    A = squeeze(DMDmode(field, :, :, :));
elseif TimePos == 5
    A = squeeze(DMDmode(field, :, :, data_shape(4), :));
end

figure(2);
nModes = length(freq);
for i = 1:nModes
    subplot(nModes, 1, i);
    contourf(real(squeeze(A(:, :, indices(i)))), 50, 'EdgeColor','none');
    title(sprintf('Freq = %.2f', freq(i)), 'FontSize', 14, 'Interpreter','latex');
end
set(gcf, "Position", [100 150 500 600])

