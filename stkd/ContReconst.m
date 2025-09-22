function ContReconst=ContReconst(t,t0,u,Tiempos,deltas,omegas)
%M=length(deltas);
[N,M]=size(u);
% size(u)
% return
% caci=zeros(N,1);
% for m=1:M
%  caci=caci+u(:,m)*exp((deltas(m)+i*omegas(m))*(t-t0));   
% end
% ContReconst=caci;

vv=zeros(M,1);
for m=1:M
 vv(m)=exp((deltas(m)+i*omegas(m))*(t-t0));   
end
%vv=vv;
ContReconst=u*vv;
%size(ContReconst)