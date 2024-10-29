function [F]=dyn(Ntotal,Ntotali,d,Mf,NN,LELE,LEPLE,F0)

F=zeros(Ntotal,Ntotal);
Fv=zeros(Ntotal^2,1);

parfor kk=1:Ntotal^2

j=Ntotali(kk,1);
i=Ntotali(kk,2);
Fv(kk)=ingF(j,i,d,Mf,NN,LELE,LEPLE,F0);

end

for kk=1:Ntotal^2
j=Ntotali(kk,1);
i=Ntotali(kk,2);
F(j,i)=Fv(kk);
end


