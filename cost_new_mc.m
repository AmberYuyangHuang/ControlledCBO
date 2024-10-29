%For Ackley function
function [L]=cost_new_mc(E)
global Ntotal d NN lim
L=zeros(Ntotal,1);

%number of sampling point for MC integration
num=10e6;
x=-lim+2*lim*rand(d,num);
for i=1:Ntotal
    sum=E(x);
    for p=1:d
        sum=sum.*leg(x(p,:),NN(i,p));
    end
L(i)=((2*lim)^d)*mean(sum);
end
