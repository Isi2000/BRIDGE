
function  [GrowthRate,Frequency,hatMode,Amplitude] =DMDd(d,hatT,Time,varepsilon1,varepsilon2)

%%%%%%%%%%%%%%%%%%%%%%%%%  DMD-d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function solves the HODMD algorithm presented in               %%%
%%% Le Clainche & Vega, SIAM J. on Appl. Dyn. Sys. 16(2):882-925, 2017  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% %% INPUT: %%
%%% d: parameter of DMD-d (higher order Koopman assumption)
%%% hatT: reduced snapshot matrix
%%% Time: vector time
%%% varepsilon1: first tolerance (SVD)
%%% varepsilon: second tolerance (DMD-d modes)
%%% %% OUTPUT: %%
%%% GrowthRate: growht rate of DMD modes
%%% Frequency: frequency of DMD modes(angular frequency)
%%% Amplitude: amplitude of DMD modes
%%% hatMode: DMD modes in reduced dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,K]=size(hatT);

%% Create the modified snapshot matrix
tildeT=zeros(d*N,K-d+1);
for ppp=1:d
    tildeT((ppp-1)*N+1:ppp*N,:)=hatT(:,ppp:ppp+K-d);
end

%% Dimension reduction
[U1,Sigma1,T1]=svd(tildeT,'econ');
sigmas1=diag(Sigma1);

Deltat=Time(2)-Time(1);
n=length(sigmas1);

NormS=norm(sigmas1,2);
kk1=0;
for k=1:n
    RRMSEE(k)=norm(sigmas1(k:n),2)/NormS;
    if RRMSEE(k)>varepsilon1
        kk1=kk1+1;
    end
end

('Spatial dimension reduction')
kk1

U1=U1(:,1:kk1);
hatT1=Sigma1(1:kk1,1:kk1)*T1(:,1:kk1)';

%% Reduced modified snapshot matrix
[~,K1]=size(hatT1);
[tildeU1,tildeSigma,tildeU2]=svd(hatT1(:,1:K1-1),'econ');

%% Reduced modified Koopman matrix
tildeR=hatT1(:,2:K1)*tildeU2*inv(tildeSigma)*tildeU1';
[tildeQ,tildeMM]=eig(tildeR);
autovalores=diag(tildeMM);

M=length(autovalores);
qq=log(autovalores);
GrowthRate=real(qq)/Deltat;
Frequency=imag(qq)/Deltat;

Q=U1*tildeQ;
Q=Q((d-1)*N+1:d*N,:);
[NN,MMM]=size(Q);

for m=1:MMM
    NormQ=Q(:,m);
    Q(:,m)= Q(:,m)/norm(NormQ(:),2);
end

%% Calculate amplitudes
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

%% Spectral complexity: number of DMD modes
('Spectral complexity')
kk3
hatMode=u(:,1:kk3);
Frequency=Frequency(1:kk3);
GrowthRate=GrowthRate(1:kk3);
Amplitude=Amplitude(1:kk3);
