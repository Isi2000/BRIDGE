%% CLEAR
clc, clear all, close all

%% LOAD DATABASE
% As a test case you may download the 2D flow past a cylinder database.
% LINK: https://drive.google.com/drive/folders/1eaL945MC46rwhsft72LE_GtSCl15ugB4

load ./../Tensor.mat % Load data tensor
load ./../X.mat      % Load mesh
load ./../Y.mat      % Load mesh

%% CALIBRATION
d_time_value = [50 100]; % Number of HODMD windows for temporal expansion
d_space_value = [1 5 10];     % Number of HODMD windows for spatial expansion
tol_value = [1e-2 1e-3];     % Tolerance to truncate the number of modes retained
dT = 1;                  % Time step size of the database
TimePos = ndims(Tensor); % Define the temporal dimension
dX = 0.07;               % Spatial grid size of the desired direction
SpacePos = 3;            % Define desired spatial direction (2, 3 or 4)

%% PERFORM MULTIDIMENSIONAL STKD
nT = size(Tensor,TimePos);
nX = size(Tensor,SpacePos);
Time = [1:nT]*dT;
Space = [1:nX]*dX;
nn0 = size(Tensor);

for tol = tol_value
    % Set tolerances
    varepsilon1=tol % SVD tol.
    varepsilon2=tol % DMD tol.
    
    % Initialize the variable containing the number of SVD modes
    nn1(1)=nn0(1);
    nn1(2:length(nn0))=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HOSVD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [hatT,U,S,sv,nn,N]=hosvd_function(Tensor,varepsilon1,TimePos);
    % hatT: Reduced temporal matrix
    % U:    Temporal SVD modes
    % S:    Tensor core
    % sv:   Singular values
    % nn:   Number of singular values
    % N:    Number of singular values in hatT (temporal SVD modes)
    
    for d_time = d_time_value
        for d_space = d_space_value
            
            %%%%%%%%%%%%%%%%%%%%%% SAVE RESULTS IN FOLDER DMD_solution %%%%%%%%%%%%%%%
            mkdir('DMD_solution')
            system('rm -r DMD_solution');
            mkdir('DMD_solution')
            filename = sprintf('./DMD_solution/DMD_history.txt' );
            diary(filename)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % DMD Temporal
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            d = d_time;
            
            if d>1
                [GrowthRateT,FrequencyT,hatModeT,AmplitudeT] =DMDdOptimReduced(d,hatT,Time,varepsilon1,varepsilon2);
            else
                [GrowthRateT,FrequencyT,hatModeT,AmplitudeT] =DMD1OptimReduced(hatT,Time,varepsilon1,varepsilon2);
            end
            
            ('HODMD complete!')
            [N,K] = size(hatT);
            hatTReconst=zeros(N,K);
            for k=1:K
                hatTReconst(:,k)= ContReconst(Time(k),Time(1),hatModeT,Time,GrowthRateT,FrequencyT);
            end
            errorT = norm(hatT-hatTReconst,2)/norm(hatT,2)
            
            GrowthRateFrequencyAmplitudeTemporal = [GrowthRateT',FrequencyT',AmplitudeT']
            save ./DMD_solution/GrowthRateFrequencyAmplitudeTemporal.mat GrowthRateFrequencyAmplitudeTemporal
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % DMD Spatial
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            d = d_space;
            
            hatX = diag(sv{SpacePos})*U{SpacePos}';
            
            if d>1
                [GrowthRateX,FrequencyX,hatModeX,AmplitudeX] =DMDdOptimReduced(d,hatX,Space,varepsilon1,varepsilon2);
            else
                [GrowthRateX,FrequencyX,hatModeX,AmplitudeX] =DMD1OptimReduced(hatX,Space,varepsilon1,varepsilon2);
            end
            
            [N,I] = size(hatX);
            
            hatXReconst=zeros(N,I);
            for k=1:I
                hatXReconst(:,k)= ContReconst(Space(k),Space(1),hatModeX,Space,GrowthRateX,FrequencyX);
            end
            errorX = norm(hatX-hatXReconst,2)/norm(hatX,2)
            
            GrowthRateFrequencyAmplitudeSpatial = [GrowthRateX',FrequencyX',AmplitudeX']
            
            save ./DMD_solution/GrowthRateFrequencyAmplitudeSpatial.mat GrowthRateFrequencyAmplitudeSpatial
            
            
            % Spatio-temporal reconstruction
            
            Unew=U;
            Umodes=U;
            HH=hatTReconst;
            PP=hatModeT;
            for kk=1:nn(TimePos)
                %UTNuevo{5}(kk,:)=hatTReconst(kk,:)/sqrt(sv{5}(kk));
                HH(kk,:)=hatTReconst(kk,:)/sv{TimePos}(kk);
                PP(kk,:)=hatModeT(kk,:)/sv{TimePos}(kk);
            end
            Unew{TimePos}=HH';
            Umodes{TimePos}=conj(PP');
            
            HH=hatXReconst;
            PP=hatModeX;
            for kk=1:nn(SpacePos)
                %UTNuevo{2}(kk,:)=hatXReconst(kk,:)/sqrt(sv{2}(kk));
                HH(kk,:)=hatXReconst(kk,:)/sv{SpacePos}(kk);
                PP(kk,:)=hatModeX(kk,:)/sv{SpacePos}(kk);
            end
            Unew{SpacePos}=HH';
            Umodes{SpacePos}=conj(PP');
            modesReconst=tprod(S,Umodes);
            
            % Reconstruct the tensor saving all spatio-temporal modes
            TensorReconst=tprod(S, Unew);
            save('./DMD_solution/TensorReconst.mat','TensorReconst','-v7.3')
            
            RRMSE_SpatioTemporal = norm(Tensor(:)-TensorReconst(:),2)/norm(Tensor(:),2)
            
            
            % Spatio-temporal amplitude matrix
            dims = size(modesReconst);
            ND = ndims(modesReconst);
            
            NN=dims(TimePos);
            MM=dims(SpacePos);
            Amplitudes = zeros(MM, NN);
            ModesNumber = 0;
            
            for ene = 1:NN  % Time
                for eme = 1:MM  % Space
                    % Create colon index for all dims
                    idx = repmat({':'}, 1, ND);
                    idx{SpacePos} = eme;
                    idx{TimePos} = ene;
            
                    aba = modesReconst(idx{:});
                    Amplitudes(eme, ene) = norm(aba(:), 2);
            
                    if Amplitudes(eme, ene) / max(Amplitudes(:)) < varepsilon2
                        Amplitudes(eme, ene) = 0;
                    else
                        ModesNumber = ModesNumber + 1;
                    end
                end
            end
            
            ('Spatial Modes')
            MM
            ('Temporal Modes')
            NN

            pepi=[Amplitudes(1:MM,1:NN)',FrequencyT(1:NN)'];
            pepi=sortrows(pepi,-(MM+1));
            ampl=pepi(:,1:MM);
            Frequencyt00=pepi(:,MM+1);
            
            pepa=[ampl',FrequencyX(1:MM)'];
            pepa=sortrows(pepa,-(NN+1));
            Amplitudes=pepa(:,1:NN)';
            Frequencyx00=pepa(:,NN+1);
            
            h5=figure;
            NumModes=0;
            for n=1:NN
                for m=1:MM
                    if Amplitudes(n,m)/max(max(Amplitudes))>varepsilon2
                        semilogy(-Frequencyx00(m)/pi,Amplitudes(n,m)/max(max(Amplitudes)),'o','color','k');
                        NumModes=NumModes+1;
                        hold on
                    end
                end
            end
            xlabel('k_m')
            ylabel('a_{mn}')
            name5 = sprintf('./DMD_solution/OmegaXAmplitudes' );
            saveas(h5,name5,'fig')
            
            
            h6=figure;
            NumModes=0;
            for n=1:NN
                for m=1:MM
                    if Amplitudes(n,m)/max(max(Amplitudes))>varepsilon2
                        semilogy(Frequencyt00(n),Amplitudes(n,m)/max(max(Amplitudes)),'o','color','k');
                        NumModes=NumModes+1;
                        hold on
                    end
                end
            end
            
            xlabel('\omega_n')
            ylabel('a_{mn}')
            name5 = sprintf('./DMD_solution/OmegaTAmplitudes' );
            saveas(h5,name5,'fig')
            
            h7=figure;
            for n=1:NN
                for m=1:MM
                    if Amplitudes(n,m)/max(max(Amplitudes))>varepsilon2
                        plot(-Frequencyx00(m)/pi,Frequencyt00(n),'o','color','k');
                        hold on
                    end
                end
            end
            
            xlabel('k_m')
            ylabel('\omega_n')
            name7 = sprintf('./DMD_solution/omegaXomegaT' );
            saveas(h7,name7,'fig')
            
            save ./DMD_solution/FrequencyX Frequencyx00
            save ./DMD_solution/FrequencyT Frequencyt00
            save ./DMD_solution/AmplitudesTX Amplitudes
            
            % RECONSTRUCTION OF DMD MODES: spatio-temporal modes
            
            idx = repmat({':'}, 1,  ND);
            idx{SpacePos} = 1:MM;
            idx{TimePos}  = 1:NN;
            
            modesReconstB = modesReconst(idx{:});
            
            % Exponentials
            ExpT = zeros(K,NN);
            ExpX = zeros(I,MM);
            for kkk = 1:K
                for ene = 1:NN
                    ExpT(kkk,ene) = exp((GrowthRateT(ene) + 1i*FrequencyT(ene)) * (Time(kkk) - Time(1)));
                end
            end
            for iii = 1:I
                for eme = 1:MM
                    ExpX(iii,eme) = exp((GrowthRateX(eme) + 1i*FrequencyX(eme)) * (Space(iii) - Space(1)));
                end
            end
            
            % Generalized Prev
            dims = nn0;
            dims(TimePos) = NN;
            Prev = zeros(dims);
            
            for iii = 1:I
                modrec_dims = nn0;
                modrec_dims(SpacePos) = 1;
                modrec_dims(TimePos) = NN;
                ModRec = zeros(modrec_dims);
            
                for eme = 1:MM
                    idx = repmat({':'}, 1, ndims(modesReconstB));
                    idx{SpacePos} = eme;
                    ModRec = ModRec + modesReconstB(idx{:}) * ExpX(iii,eme);
                end
            
                idx_prev = repmat({':'}, 1, ndims(Prev));
                idx_prev{SpacePos} = iii;
                Prev(idx_prev{:}) = ModRec;
            end
            
            % TensorReconstReduced
            dimsTR = size(Prev);
            dimsTR(TimePos) = K;
            TensorReconstReduced = zeros(dimsTR);
            
            for kkk = 1:K
                modrec_dims = size(Prev);
                modrec_dims(TimePos) = 1;
                ModRec = zeros(modrec_dims);
            
                for ene = 1:NN
                    idx_prev = repmat({':'}, 1, ndims(Prev));
                    idx_prev{TimePos} = ene;
                    ModRec = ModRec + Prev(idx_prev{:}) * ExpT(kkk,ene);
                end
            
                idx_tr = repmat({':'}, 1, ndims(TensorReconstReduced));
                idx_tr{TimePos} = kkk;
                TensorReconstReduced(idx_tr{:}) = ModRec;
            end
            
            % Error computation
            errorTensorReduced = norm(TensorReconstReduced(:) - Tensor(:),2) / norm(Tensor(:),2);
            disp(['Relative error: ', num2str(errorTensorReduced)])
            
            save('./DMD_solution/TensorReconstReduced.mat','TensorReconstReduced','-v7.3');
            
            % DMD mode
            DMDmode = modesReconstB;
            save('./DMD_solution/DMDmode_tensor.mat', 'DMDmode');
            
            % DMD mode temporal (collapse spatial modes)
            dimsTemp = size(modesReconstB);
            dimsTemp(SpacePos) = I;
            DMDmode_Temporal = zeros(dimsTemp);
            
            for iii = 1:I
                modrec_dims = size(modesReconstB);
                modrec_dims(SpacePos) = 1;
                ModRec = zeros(modrec_dims);
            
                for eme = 1:MM
                    idx = repmat({':'}, 1, ndims(modesReconstB));
                    idx{SpacePos} = eme;
                    ModRec = ModRec + modesReconstB(idx{:}) * ExpX(iii,eme);
                end
            
                idx_temp = repmat({':'}, 1, ndims(DMDmode_Temporal));
                idx_temp{SpacePos} = iii;
                DMDmode_Temporal(idx_temp{:}) = ModRec;
            end
            save('./DMD_solution/DMDmode_Temporal.mat', 'DMDmode_Temporal');
            
            
            % DMD mode spatial (collapse temporal modes)
            dimsSpat = size(modesReconstB);
            dimsSpat(TimePos) = K;
            DMDmode_Spatial = zeros(dimsSpat);
            
            for iii = 1:K
                modrec_dims = size(modesReconstB);
                modrec_dims(TimePos) = 1;
                ModRec = zeros(modrec_dims);
            
                for ene = 1:NN
                    idx = repmat({':'}, 1, ndims(modesReconstB));
                    idx{TimePos} = ene;
                    ModRec = ModRec + modesReconstB(idx{:}) * ExpT(iii,ene);
                end
            
                idx_spat = repmat({':'}, 1, ndims(DMDmode_Spatial));
                idx_spat{TimePos} = iii;
                DMDmode_Spatial(idx_spat{:}) = ModRec;
            end
            save('./DMD_solution/DMDmode_Spatial.mat', 'DMDmode_Spatial');
            
            
            diary off
            close all
            
            name_folder=sprintf('Cases/mdSTKD_solution_dTime%d_dSpace%d_tol%0.1e',d_time,d_space,tol)
            movefile('DMD_solution',name_folder)
        end
    end
end