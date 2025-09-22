
function  [deltas,omegas,hatmodos,hatamplitudes] =DMDdOptimReduced(d,hatT,Tiempos,varepsilon1,varepsilon2)
                
[Jprima,K]=size(hatT);

N=Jprima;

tildeT=zeros(d*N,K-d+1);
%size(tildeT)
for ppp=1:d
 tildeT((ppp-1)*N+1:ppp*N,:)=hatT(:,ppp:ppp+K-d);   
end
% size(tildeT)
% return
[U1,Sigma1,T1]=svd(tildeT);
sigmas1=diag(Sigma1);
%sigmas
Deltat=Tiempos(2)-Tiempos(1);
n=length(sigmas1);

aaaaa=norm(sigmas1,2);
kk1=0;
for k=1:n
    RRMSEE(k)=norm(sigmas1(k:n),2)/aaaaa;
        if RRMSEE(k)>varepsilon1
        kk1=kk1+1;
        end
end

n
kk1
%kk1=8
close(figure(2))
figure(2)
%semilogy(1:kk1,sigmas1(1:kk1));
semilogy(1:n,RRMSEE(1:n)/RRMSEE(1),'k','linewidth',1);
hold on
semilogy(1:kk1,RRMSEE(1:kk1)/RRMSEE(1),'k','linewidth',3);
%close(figure(1))

U1=U1(:,1:kk1);

hatT1=Sigma1(1:kk1,1:kk1)*T1(:,1:kk1)';
% K
% return
[~,K1]=size(hatT1);


[tildeU1,tildeSigma,tildeU2]=svd(hatT1(:,1:K1-1),'econ');

tildeR=hatT1(:,2:K1)*tildeU2*inv(tildeSigma)*tildeU1';
[tildeQ,tildeMM]=eig(tildeR);
autovalores=diag(tildeMM);



M=length(autovalores);
cucu=log(autovalores);
deltas=real(cucu)/Deltat;
omegas=imag(cucu)/Deltat;
%omegas
Q=U1*tildeQ;

Q=Q((d-1)*N+1:d*N,:);

%Q=Q(1:N,:);

 [NN,MMM]=size(Q);
 
 for m=1:MMM
    aaaaa=Q(:,m);
   Q(:,m)= Q(:,m)/norm(aaaaa(:),2);
 end

meme=zeros(NN*K,M);
bebe=zeros(NN*K,1);
aa=eye(MMM);
for k=1:K
 meme(1+(k-1)*NN:k*NN,:)=Q*aa; 
 aa=aa*tildeMM;
 bebe(1+(k-1)*NN:k*NN,1)=hatT(:,k);
end

[Ur,Sigmar,Vr]=svd(meme,'econ');
a=Vr*(Sigmar\(Ur'*bebe));
%return

u=zeros(NN,M);
for m=1:M
    u(:,m)=a(m)*Q(:,m);
end
hatamplitudes=zeros(M,1);

for m=1:M
    hh=deltas(m);
    pepi=exp(hh.*Tiempos);
    pepe=norm(pepi,2)/sqrt(K);
    hatamplitudes(m)=norm(u(:,m),2)*pepe;    
end


UU=[u;deltas';omegas';hatamplitudes']';
UU1=sortrows(UU,-(NN+3));

UU=UU1';
u=UU(1:NN,:);
deltas=UU(NN+1,:);
omegas=UU(NN+2,:);
hatamplitudes=UU(NN+3,:);
% hatamplitudes
% return
kk3=0;
for m=1:M
    if hatamplitudes(m)/hatamplitudes(1)>varepsilon2
        kk3=kk3+1;
    else
    end
end
kk3

%amplitudes=UU(NN+3,:);
hatmodos=u(:,1:kk3);
omegas=omegas(1:kk3);
deltas=deltas(1:kk3);
hatamplitudes=hatamplitudes(1:kk3);
