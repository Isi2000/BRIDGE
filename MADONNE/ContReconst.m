function Reconst=ContReconst(t,t0,u,Time,GrowthRate,Frequency)

[N,M]=size(u);
vv=zeros(M,1);
for m=1:M
 vv(m)=exp((GrowthRate(m)+i*Frequency(m))*(t-t0));   
end

Reconst=u*vv;
