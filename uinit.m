function [U]=uinit(Nu)
global Ntotal NN gamma d LE gp gw
U=zeros(Ntotal,1);
if (Nu>0)
for i=1:Ntotal
    wa=0;   
 for l=1:Nu
    for j=1:Nu
        prod=1;
        for k=1:d
            %LE:phi_i^k(x_k)  NN:basis permutation
            fact=gw*(u0sep(l,k,gp).*u0sep(j,k,gp).*LE(:,NN(i,k)+1));
            %product over k
            prod=prod*fact;
        end
        wa=wa+prod;
    end
 end
 U(i)=gamma.*wa;
end
end