%% CLEAR
clc, clear all, close all
%% LOAD DATABASE
% As a test case you may download the 2D flow past a cylinder database.
% LINK: https://drive.google.com/drive/folders/1eaL945MC46rwhsft72LE_GtSCl15ugB4

load /home/isacco/DATASETS/BINNED_TEMPERATURES/binned_temperature_cs.mat % Load data tensor
start_t = 1;
end_t = 100;
Tensor = binned_temperature_cs(:, :, :, start_t:end_t); % Assign the loaded data to the Tensor variable
clear binned_temperature_cs
%% CALIBRATION
d_value = [10 25 40 ];      % Number of HODMD windows (Recommendation 0.1-0.4*Nt)
tol_value = [1e-3 1e-4];   % Tolerance to truncate the number of modes retained
dT = 1;                % Time step of the database
n_iter = 1;            % Change to 1 if non-iterative approach is wanted

%% PERFORM MULTIDIMENSIONAL ITERATIVE HODMD
TimePos = ndims(Tensor);
Nt = size(Tensor,TimePos);
Time = [1:Nt]*dT;
nn0 = size(Tensor);

for tol = tol_value
    
    clear hat* DMD* Mode* ome* del* ampl
    
    varepsilon1=tol % Tolerance for SVD
    varepsilon2=tol % Tolerance for DMD

    % Initialize the variable containing the number of SVD modes
    nn1(1)=nn0(1);
    nn1(2:length(nn0))=0;

    % PERFORM HOSVD
    [hatT,U,S,sv,nn,N]=hosvd_function(Tensor,varepsilon1,TimePos);
    % hatT: Reduced temporal matrix
    % U:    Temporal SVD modes
    % S:    Tensor core
    % sv:   Singular values
    % nn:   Number of singular values
    % N:    Number of singular values in hatT (temporal SVD modes)

    for d = d_value
        
        mkdir('DMD_solution')
        system('rm -r DMD_solution');
        mkdir('DMD_solution')
        filename = sprintf('./DMD_solution/DMD_history.txt' );
        diary(filename)

        % PERFORM DMD-d TO THE REDUCED TEMPORAL MATRIX hatT  
        if d>1
            [GrowthRate,Frequency,hatMode,Amplitude] = DMDd(d,hatT,Time,varepsilon1,varepsilon2);
        else
            [GrowthRate,Frequency,hatMode,Amplitude] = DMD1(hatT,Time,varepsilon1,varepsilon2);
        end
        % GrowthRate:   Growth rate of the DMD modes
        % Frequency:    Frequency of the DMD modes
        % hatMode:      Reduced DMD mode
        % Amplitude:    Amplitudes weighting the DMD modes        

        % Reconstruct the original Tensor using the DMD expansion:
        Tensor_Reconst = DMDreconst(GrowthRate,Frequency,hatMode,Time,U,S,sv,hatT,nn);

        % ITERATIVE LOOP
        for zz = 2:n_iter

		    ('Iteration number:' + string(zz))
		
		    TensorIT = Tensor_Reconst;

            ('Performing HOSVD. Please wait...')
		    [hatTIT,UIT,SIT,svIT,nn1IT,NIT]=hosvd_function(Tensor,varepsilon1,TimePos);
            ('HOSVD complete!  ')
		
		    
		    ('Performing HODMD. Please wait...')
		    if d>1
                [GrowthRate,Frequency,hatMode,Amplitude] = DMDd(d,hatTIT,Time,varepsilon1,varepsilon2);
            else
                [GrowthRate,Frequency,hatMode,Amplitude] = DMD1(hatTIT,Time,varepsilon1,varepsilon2);
            end
		    ('HODMD complete!')
		    
		    Tensor_Reconst = DMDreconst(GrowthRate,Frequency,hatMode,Time,UIT,SIT,svIT,hatTIT,nn1IT);
		    
		    % Break the loop when the number of singular values is the same in two consecutive iterations:
		    if nn1IT==nn
			    break
            end
		    nn = nn1IT;

        end
        
        if n_iter>1
            sprintf('Final number of iterations for d = %d and tol = %1.1i: %d',d,tol,zz-1)
        end

        % Don't do this for complex data
        Tensor_Reconst = real(Tensor_Reconst); % Reconstruct database
        
        ('Relative mean square error made in the calculations')
        RRMSE = norm(Tensor(:)-Tensor_Reconst(:),2)/norm(Tensor(:),2)

        ('Growth rate, Frequency, Amplitude')
        GrowthRateFrequencyAmplitude = [GrowthRate',Frequency',Amplitude']
               
        diary off
        
        % Save reconstruction and mode information
        save ./DMD_solution/Tensor_Reconst.mat Tensor_Reconst -v7.3
        save ./DMD_solution/GrowthrateFrequencyAmplitude.mat GrowthRateFrequencyAmplitude
        

        % CALCULATE DMD MODES
        ('Calculating DMD modes...')
        [N,~]=size(hatT);
        hatMode_m=zeros(N,length(Amplitude));
        for ii=1:length(Amplitude)
            hatMode_m(:,ii)=hatMode(:,ii)/Amplitude(ii);
        end
        
        % Temporal DMD modes in reduced dimension
        Modes = U;
        Modes{TimePos} = hatMode_m';
        % Reconstruction of the temporal DMD modes
        DMDmode=tprod(S,Modes);
        % Save DMD modes
        save ./DMD_solution/DMDmode.mat DMDmode -v7.3
        
        name_folder=sprintf('Cases/mdHODMDit_solution_d%d_tol%0.1e',d,tol)
        movefile('DMD_solution',name_folder)

    end
end
