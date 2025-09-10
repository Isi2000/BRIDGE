function DMDmode=calculateDMDmode(N,hatMode,Amplitude,U,S,nn,TimePos)

hatMode_m=zeros(N,length(Amplitude));
for ii=1:length(Amplitude)
    hatMode_m(:,ii)=hatMode(:,ii)/Amplitude(ii);
end

for kk=1:nn(TimePos)
    ModesT(kk,:)=hatMode_m(kk,:);
end
% Temporal DMD modes in reduced dimension
Modes=U;
Modes{TimePos}=ModesT';
% Reconstruction of the temporal DMD modes
DMDmode=tprod(S,Modes);
