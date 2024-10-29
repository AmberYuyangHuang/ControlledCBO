function [val]=bb()
global LELE Ntotal d NN
b=zeros(Ntotal,1);

for i=1:Ntotal
    fac=1;
    for k=1:d
        fac=fac*sum(LELE(:,NN(i,k)+1,NN(i,k)+1));
    end
    b(i)=fac;
    
end
val=b;