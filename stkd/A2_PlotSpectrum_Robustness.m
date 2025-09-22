%% CLEAR
clc, clear all, close all

%% CALIBRATION
d_time_value = [50 100];    % Number of HODMD windows in time
d_space_value = [1];      % Number of HODMD windows in space
tol_value = [1e-3];    % Tolerance to truncate the number of modes retained

%% PLOT SPECTRUM
colors = hsv(length(d_time_value)*length(d_space_value));
mark = ['o' 'x' '+' '^' '*' 's'];
j = 0; i = 0;

figure(1); hold on
tiledlayout(1,2,'TileSpacing','compact')

for tol = tol_value
    
    j = j+1; i = 0;
    
    for d_time = d_time_value
        for d_space = d_space_value
        
            i = i+1;
            
            name_folder=sprintf('./Cases/mdSTKD_solution_dTime%d_dSpace%d_tol%0.1e/GrowthrateFrequencyAmplitudeTemporal.mat',d_time,d_space,tol);
            load(name_folder)
            CT= GrowthRateFrequencyAmplitudeTemporal;
            omegasT = CT(:,2);
            amplitudesT = CT(:,3);

            name_folder=sprintf('./Cases/mdSTKD_solution_dTime%d_dSpace%d_tol%0.1e/GrowthrateFrequencyAmplitudeSpatial.mat',d_time,d_space,tol);
            load(name_folder)
            CX = GrowthRateFrequencyAmplitudeSpatial;
            omegasX = CX(:,2);
            amplitudesX = CX(:,3);
    
            name_legend = sprintf('d_t = %0.0i, d_x = %0.0i, tol  = %0.1e',d_time,d_space,tol);

            nexttile(1); hold on
            plot(omegasT,amplitudesT/max(amplitudesT),'LineStyle','none','Marker',mark(j),'Color',colors(i,:),'DisplayName',name_legend)
            set(gca,'YScale','log')
            % xlim([0 1.6]); ylim([1e-1 1])
            xlabel('\omega_m'); ylabel('a_m')
            grid  minor

            nexttile(2); hold on
            plot(omegasX,amplitudesX/max(amplitudesX),'LineStyle','none','Marker',mark(j),'Color',colors(i,:),'DisplayName',name_legend)
            set(gca,'YScale','log')
            % xlim([0 1.6]); ylim([1e-1 1])
            xlabel('\omega_m'); ylabel('a_m')
            grid  minor
            legend('Location','northeastoutside')

        end
    end 
end

set(gcf,'Position',[0 100 2000 400])
