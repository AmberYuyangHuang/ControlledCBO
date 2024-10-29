function [G]=grad0(Ntotali,Ntotal,NN,Nu,LELE,LEPLE,gp,d)
%G matrix of size Ntotal*Ntotal
N=Ntotal;
G=zeros(N,N);
Gv=zeros(N^2,1);

parfor kk=1:N^2
j=Ntotali(kk,1);
i=Ntotali(kk,2);
Gv(kk)=inG(i,j,Nu,LELE,LEPLE,gp,d,NN);
end

for kk=1:N^2
j=Ntotali(kk,1);
i=Ntotali(kk,2);
G(j,i)=Gv(kk);
end
