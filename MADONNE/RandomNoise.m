clc
clear all
close all


load ./../reduced_domain_matrix_Re220p6_t2900_3d/Tensor_Re220p6_t2900_3d.mat

%% Percentage noise level
NoiseLevel=10

%%
MaxUz=0.15
[L I J K T]=size(Tensor);


for l=1:L
    for i=1:I
        for j=1:J
            for k=1:K
                for t=1:T
                    TensorNoise(l,i,j,k,t)=Tensor(l,i,j,k,t)+rand(1)*NoiseLevel/100;
                end
            end
        end
    end
end

save Tensor_Re220p6_t2900_3d_noise.mat TensorNoise

figure(1)
contourf(squeeze(TensorNoise(1,:,:,1,10)))
figure(2)
contourf(squeeze(TensorNoise(3,:,:,1,10)))
