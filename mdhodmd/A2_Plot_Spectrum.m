%% CLEAR
clc; clear all; close all;

%% CALIBRATION
d_value = [10 25 40]; % Number of HODMD windows
tol_value = [1e-3 1e-4];      % Tolerance to truncate the number of modes retained

%% PLOT SPECTRUM
colors = ['m' 'c' 'r' 'g' 'b' 'k'];
mark = ['o' 'x' '+'];
j = 0; i = 0;

figure(1); hold on

for tol = tol_value
    
    j = j+1; i = 0;
    
    for d = d_value
        
        i = i+1;
        
        name_folder=sprintf('./Cases/mdHODMDit_solution_d%d_tol%0.1e/GrowthrateFrequencyAmplitude.mat',d,tol);
        load(name_folder)
        C1 = GrowthRateFrequencyAmplitude;
        omegas1 = C1(:,2);
        amplitudes1 = C1(:,3);

        name_legend = sprintf('d = %0.0i, tol  = %0.1e',d,tol);
        plot(omegas1,amplitudes1/max(amplitudes1),'LineStyle','none','Marker',mark(j),'Color',colors(i),'DisplayName',name_legend)

    end 
end

legend
set(gca,'YScale','log')
% xlim([0 1.6]); ylim([1e-1 1])
xlabel('\omega'); ylabel('amplitud')
grid  minor

set(gcf,'Position',[0 100 2000 400])
set(legend,'Position',[0.89 0.48 0.1 0.1])
set(gca,'Position',[0.05 0.15 0.82 0.8])