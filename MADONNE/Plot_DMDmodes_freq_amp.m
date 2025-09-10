clc
clear all
close all

tolit=[1e-5 1e-4 ]
ddit=[30 40 50 60]


figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
box on
set(axes1,'YMinorTick','on','YScale','log');
% xlim([-0.01 0.7])%
for i=1:length(tolit)
    
    for j=1:length(ddit)
        
        d=ddit(j)
        varepsilon1=tolit(i)
        a=sprintf('DMD_solution_d%0.0i_tol%0.0e/dataDeltasOmegasAmplTemporal.mat',d,varepsilon1)
        m=load(a);
        if i==1
            plot(m.DeltasOmegAmplTemporal(:,3)/DeltaTparam,m.DeltasOmegAmplTemporal(:,4)/(max(m.DeltasOmegAmplTemporal(:,4))),'s')
        elseif i==2
            plot(m.DeltasOmegAmplTemporal(:,3)/DeltaTparam,m.DeltasOmegAmplTemporal(:,4)/(max(m.DeltasOmegAmplTemporal(:,4))),'^')
        else
            plot(m.DeltasOmegAmplTemporal(:,3)/DeltaTparam,m.DeltasOmegAmplTemporal(:,4)/(max(m.DeltasOmegAmplTemporal(:,4))),'o')
        end
    end
end


xlabel('\omega_m')
ylabel('a_{m}')
%ylim([1e-3 1])
%xlim([-0.01 6])
%legend('d12-Inf','d12-L2','d15-Inf','d15-L2','d18-Inf','d18-L2','d20-Inf','d20-L2')
set(axes1,'FontSize',18,'YMinorTick','on','YScale','log');
set(gcf, 'Position', [100, 100, 1000, 300])
 
