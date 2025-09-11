
function  [GrowthRate,Frequency,hatMode,Amplitude] =DMD1(hatT,Time,varepsilon1,varepsilon2)

%%%%%%%%%%%%%%%%%%%%%%%%%  DMD-1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function solves the DMD-1 algorithm presented in               %%%
%%% Le Clainche & Vega, SIAM J. on Appl. Dyn. Sys. 16(2):882-925, 2017  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% %% INPUT: %%
%%% hatT: reduced snapshot matrix
%%% Time: vector time
%%% varepsilon1: first tolerance (SVD)
%%% varepsilon: second tolerance (DMD-d modes)
%%% %% OUTPUT: %%
%%% GrowthRate: growht rate of DMD modes
%%% Frequency: frequency of DMD modes(angular frequency)
%%% Amplitude: amplitude of DMD modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Create reduced snapshot matrix
[N,K]=size(hatT);
[hatU1,hatSigma,hatU2]=svd(hatT(:,1:K-1),'econ');
sigmas=diag(hatSigma);

Deltat=Time(2)-Time(1);

NormS=norm(sigmas,2);
kk1=0;
n=length(sigmas);
for k=1:n
    RRMSEE=norm(sigmas(k:n),2)/NormS;
        if RRMSEE>varepsilon1
        kk1=kk1+1;
        end
end
('Spatial complexity')
kk1
hatU1=hatU1(:,1:kk1);
hatU2=hatU2(:,1:kk1);
hatSigma=hatSigma(1:kk1,1:kk1);

%% Calculate Koopman operator
hatR=hatT(:,2:K)*hatU2*inv(hatSigma)*hatU1';
[Q,tildeMM]=eig(hatR);
eigenvalues=diag(tildeMM);

M=length(eigenvalues);
qq=log(eigenvalues);
GrowthRate=real(qq)/Deltat;
Frequency=imag(qq)/Deltat;

%% Calculate amplitudes
 [NN,MMM]=size(Q);
 
 for m=1:MMM
    EigN=Q(:,m);
   Q(:,m)= Q(:,m)/norm(EigN(:),2);
 end

Mm=zeros(NN*K,M);
Bb=zeros(NN*K,1);
aa=eye(MMM);
for k=1:K
 Mm(1+(k-1)*NN:k*NN,:)=Q*aa; 
 aa=aa*tildeMM;
 Bb(1+(k-1)*NN:k*NN,1)=hatT(:,k);
end

[Ur,Sigmar,Vr]=svd(Mm,'econ');
a=Vr*(Sigmar\(Ur'*Bb));
%return

u=zeros(NN,M);
for m=1:M
    u(:,m)=a(m)*Q(:,m);
end
Amplitude=zeros(M,1);

for m=1:M
    GR=GrowthRate(m);
    AmplGR=exp(GR.*Time);
    AmplN=norm(AmplGR,2)/sqrt(K);
    Amplitude(m)=norm(u(:,m),2)*AmplN;    
end

UU=[u;GrowthRate';Frequency';Amplitude']';
UU1=sortrows(UU,-(NN+3));

UU=UU1';
u=UU(1:NN,:);
GrowthRate=UU(NN+1,:);
Frequency=UU(NN+2,:);
Amplitude=UU(NN+3,:);

%% Set the number of DMD modes
kk3=0;
for m=1:M
    if Amplitude(m)/Amplitude(1)>varepsilon2
        kk3=kk3+1;
    else
    end
end
('Spectral complexity')
kk3
hatMode=u(:,1:kk3);
Frequency=Frequency(1:kk3);
GrowthRate=GrowthRate(1:kk3);
Amplitude=Amplitude(1:kk3);
