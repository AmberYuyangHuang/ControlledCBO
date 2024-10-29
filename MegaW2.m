function [val]=MegaW2(j,k,i,LELELE,LEPLEPLE)
global NN d
U=zeros(d,1);
for l=1:d
    prod=1;
for m=1:d
   switch m
        case l
        fact=LEPLEPLE(NN(j,l)+1,NN(k,l)+1,NN(i,l)+1);
        otherwise
        fact=LELELE(NN(j,m)+1,NN(k,m)+1,NN(i,m)+1);
    end
    prod=prod*fact;
end
U(l)=prod;
end

val=sum(U);

%