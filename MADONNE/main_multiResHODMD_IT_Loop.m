
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%         This Matlab file uses the Multi-dimensional HODMD           %%%
%%%      Extension of HODMD for the analysis of complex signals,        %%%
%%%                   noisy data, complex flows...                      %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The algorithm and one application for complex flows is presented in:   %
%                                                                         %
%    Le Clainche, S., Vega, J.M. & Soria, J., Higher order dynamic mode   %
%    decomposition of noisy experimental data: The flow structure of a    %
%   zero-net-mass-flux jet,Exp. Therm. Fluid Sci., 2018, 88, pp. 336-353  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   HODMD by Le Clainche, S. & Vega, J.M., Higher order dynamic mode  %%%
%%%      decomposition, SIAM J. Appl. Dyn. Sys., 16(2), pp. 882-925     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  This code presents an example describing the near field of the wake of %
%                a three-dimensional cylinder at Re=220                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% %% INPUT: %%
%%% d: parameter of DMD-d (higher order Koopman assumption)
%%% V: snapshot matrix
%%% deltaT: time between snapshots
%%% varepsilon1: first tolerance (SVD)
%%% varepsilon: second tolerance (DMD-d modes)
%%% %% OUTPUT: %%
%%% Vreconst: reconstruction of the snapshot matrix V
%%% deltas: growht rate of DMD modes
%%% omegas: frequency of DMD modes(angular frequency)
%%% amplitude: amplitude of DMD modes
%%% modes: DMD modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear all
close all
%%%%%%%%%%%%%%%%%%%%%% SAVE RESULTS IN FOLDER DMD_solution %%%%%%%%%%%%%%%
mkdir('DMD_solution')
system('rm -r DMD_solution');
mkdir('DMD_solution')
filename = sprintf('./DMD_solution/DMD_history.txt' );
diary(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% USER:
% load data in tensor form: Tensor_lij(s)k:
% l: components (i.e.:streamwise, normal and spanwise velocity)
% i: X (or Y)
% j: Y (or X)
% s: Z
% k: Time
load ./../reduced_domain_matrix_Re220p6_t2900_3d/Tensor_Re220p6_t2900_3d.mat

%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the snapshot number
SNAP=170
%varepsilon1=1e-2
%varepsilon2=1e-2
%d=1
deltaT=1;
% Set the position of the temporal variable. The code is prepared to set
% the temporal evolution in the last component of the tensor, so this
% variable also determines the dimension of the tensor
TimePos=4%5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Tensor dimension - number of snapshots:
Tensor0=Tensor;
clear Tensor
if TimePos==5
    Tensor(:,:,:,:,:)=Tensor0(:,:,:,:,1:SNAP);
elseif TimePos==4
    Tensor(:,:,:,:)=Tensor0(:,:,:,1:SNAP);
end
Time=[1:SNAP]*deltaT;
%% ALGORITHM:
%% ITERATIVE
nn0=size(Tensor);
nn(1)=nn0(1);
nn(2:length(nn0))=0;

% INTRODUCE tolerance: tolit=varepsilon1=varepsilon2 & d=ddit
tolit=[1e-4 1e-3 1e-2]
ddit=[10 20 30 40 50 60 70]
for ii=1:length(tolit)
    for jj=1:length(ddit)
        d=ddit(jj)
        varepsilon1=tolit(ii)
        varepsilon2=varepsilon1
        %
        for zz=1:1000
            ('Iteration number')
            zz
            if zz~=1
                clear S* U* Frequency GrowthRate hat* sv*
                %load ./DMD_solution/TensorReconst.mat
                %       load ./DMD_solution_tensor_time/dataTimeClean TimeClean
                Tensor=TensorReconst;
                %        Time=TimeClean;
            end
            
            %% PERFORM HOSVD DECOMPOSITION TO CALCULATE THE REDUCED TEMPORAL MATRIX hatT
            [hatT,U,S,sv,nn1]=hosvd_function(Tensor,varepsilon1,nn,nn0,TimePos);
            % hatT: reduced temporal matrix
            % U: temporal SVD modes
            % S: tensor core
            % sv: singular values
            % nn: number of singular values
            %
            % Tensor: initial data
            % n: dimension of Tensor
            % varepsilon1: tolerance SVD
            
            %% PERFORM DMD-d TO THE REDUCED TEMPORAL MATRIX hatT
            if d>1
                [GrowthRate,Frequency,Amplitude,hatMode] =DMDd(d,hatT,Time,varepsilon1,varepsilon2);
            else
                [GrowthRate,Frequency,Amplitude,hatMode] =DMD1(hatT,Time,varepsilon1,varepsilon2);
            end
            %
            % GrowthRate: growth rate of the DMD modes
            % Frequency: frequency of the DMD modes
            % hatMode: reduced DMD mode
            % Amplitude: amplitudes weighting the DMD modes
            % varepsilon1: tolerance for the dimension reduction via SVD
            % varepsilon2: tolerance to set the DMD modes (the amplitude of the DMD
            % modes retained is > varepsilon2
            
            %% RECONSTRUCT THE ORIGINAL TENSOR USING THE DMD EXPANSION
            [TensorReconst]=DMDreconst(GrowthRate,Frequency,hatMode,Time,U,S,sv,nn1,TimePos);
            %
            % TensorReconst: reconstruction of the original tensor
            % In the analysis of complex data comment the following line:
            TensorReconst=real(TensorReconst);
            
            %
            ('Relative mean square error made in the calculations')
            RRMSE=norm(Tensor(:)-TensorReconst(:),2)/norm(Tensor(:),2);
            %
            ('Growth rate, Frequency, Amplitude')
            GrowthrateFrequencyAmplitude=[GrowthRate',Frequency',Amplitude']
            
            % Break the loop when the number of singular values is the same in two consecutive iterations
            num=0
            for i=2:length(nn1)
                if nn(i)==nn1(i)
                    num=num+1;
                end
            end
            if num==length(nn1)-1
                break
            end
            nn=nn1;
            
        end
        
        ('Final number of iterations')
        zz
        diary off
        
        % Save the reconstruction of the tensor and the Growth rates, frequencies
        % and amplitudes
        save ./DMD_solution/TensorReconst.mat TensorReconst -v7.3
        save ./DMD_solution/GrowthrateFrequencyAmplitude.mat GrowthrateFrequencyAmplitude
        
        ('Calculating DMD modes...')
        %% Calculate DMD modes
        [N,~]=size(hatT);
        [DMDmode]=calculateDMDmode(N,hatMode,Amplitude,U,S,nn1,TimePos);
        % Save DMD modes
        save ./DMD_solution/DMDmode.mat DMDmode -v7.3
        
        %% Plot the results
        % Frequency vs. absolute value of Growth rate
        h=figure;
        semilogy(Frequency,abs(GrowthRate),'o','linewidth',2,'color','k');
        name1 = sprintf('./DMD_solution/FrequencyGrowthrate_d%03i',d );
        xlabel('\omega_m')
        ylabel('|\delta_m|')
        saveas(h,name1,'fig')
        
        % Frequency vs. amplitudes
        h2=figure;
        semilogy(Frequency,Amplitude,'o','linewidth',2,'color','k');
        name2 = sprintf('./DMD_solution/FrequencyAmplitude_d%03i',d );
        xlabel('\omega_m')
        ylabel('a_m')
        saveas(h2,name2,'fig')
        
        ('Check the folder DMD_solution!')
        a=sprintf('DMD_solution_d%0.0i_tol%0.0e',d,varepsilon1)
        movefile('DMD_solution',a)
        
    end
end

    

