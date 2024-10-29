function [We]=Weg(Megawind,Megadim,LELELE,LEPLEPLE)
global Ntotal
We=zeros(Ntotal,Ntotal,Ntotal);
Wtv=ones(Megadim,1);

M1=Megawind(:,1);
M2=Megawind(:,2);
M3=Megawind(:,3);

for kk=1:Megadim
i=M1(kk,1);
k=M2(kk,1);
j=M3(kk,1);
Wtv(kk,1)=MegaW2(j,k,i,LELELE,LEPLEPLE);
end


for kk=1:Megadim
i=M1(kk);
k=M2(kk);
j=M3(kk);

We(j,k,i)=Wtv(kk,1);
end


