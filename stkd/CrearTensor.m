
function  CrearTensor

%Este programa aplica espacio-temporal; variables mal escaladas porque phi
%carga las componentes phi, psi y T generadas con read_snapshots_sept y
%las organiza en un tensor
%Los nombres del tensor generado son

load('componentephi1.mat','phi1')
load('componentepsi1.mat','psi1')
load('componenteT1.mat','T1')
[M,I,J,K]=size(T1);
Tensor1=zeros(M,I,J,3,K);
 Tensor1(:,:,:,1,:)=phi1;
 Tensor1(:,:,:,2,:)=psi1;
 Tensor1(:,:,:,3,:)=T1;
 
 save('Tensor1.mat','Tensor1','-v7.3')
