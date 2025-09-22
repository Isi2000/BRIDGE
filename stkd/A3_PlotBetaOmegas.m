%% CLEAR
clc, clear all, close all

%% CALIBRATION
d_time_value = [50];    % Number of HODMD windows in time
d_space_value = [1];      % Number of HODMD windows in space
tol_value = [1e-3];    % Tolerance to truncate the number of modes retained

%% LOAD CASE
name_folder=sprintf('./Cases/mdSTKD_solution_dTime%d_dSpace%d_tol%0.1e/FrequencyT.mat',d_time_value,d_space_value,tol_value);
load(name_folder);
name_folder=sprintf('./Cases/mdSTKD_solution_dTime%d_dSpace%d_tol%0.1e/FrequencyX.mat',d_time_value,d_space_value,tol_value);
load(name_folder);
name_folder=sprintf('./Cases/mdSTKD_solution_dTime%d_dSpace%d_tol%0.1e/AmplitudesTX.mat',d_time_value,d_space_value,tol_value);
load(name_folder);

%% PLOT SPECTRUM TIME AND SPACE
[X,T] = meshgrid(Frequencyx00,Frequencyt00);
Amplitudes(Amplitudes==0) = NaN;

figure(1); hold on
% scatter(X(:), T(:), 36, Amplitudes(:), "filled")
scatter3(X(:), T(:), Amplitudes(:),36, Amplitudes(:), "filled")
colormap("turbo")
c = colorbar;
set(gca,'ColorScale','log')