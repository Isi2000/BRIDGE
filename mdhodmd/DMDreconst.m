function [TensorReconst]=...
    DMDreconst(GrowthRate,Frequency,hatMode,Time,U,S,sv,hatT,nn)

% Initialize variables
TimePos = length(nn);
[N,K]=size(hatT);
hatTReconst=zeros(N,K);

% Reconstruction using the DMD expansion
for k=1:K
    hatTReconst(:,k)= ContReconst(Time(k),Time(1),hatMode,Time,GrowthRate,Frequency);
end

% Reconstruction of the original tensor using the reduced tensor and the
% tensor core
Unondim=U;
for kk=1:nn(TimePos)
    UTnondim{TimePos}(kk,:)=hatTReconst(kk,:)/sv{TimePos}(kk);
end
Unondim{TimePos}=UTnondim{TimePos}';

TensorReconst=tprod(S, Unondim);

