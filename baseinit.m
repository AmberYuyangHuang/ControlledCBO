function [Ntotal,NN,Ntotali]=baseinit(N,gradmax,d)

%list all permutation of polynomial basis
wa=polycomb(0:N(1));
for i=2:d
    wa=polycomb(wa,0:N(i));
end
NN=wa';
NN(1,:)=[];
Ntotal=size(NN,1);

%compute total degree
grados=sum(NN,2);
%find the basis statify the max-grad requirement
lista=[];

for j=1:Ntotal
    if (grados(j)<=gradmax) 
        lista=[lista j];
    end
end
%the target basis
NN2=NN(lista,:);
NN=NN2;
%total number of possible basis
Ntotal=length(NN);

Ntotali=polycomb(1:Ntotal,1:Ntotal)';
