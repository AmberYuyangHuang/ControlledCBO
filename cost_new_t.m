
function [L]=cost_new_t()
global Ntotal d LE NN F0 Nf gw

L=zeros(Ntotal,1);
for i=1:Ntotal
suma=0;
for l=1:Nf
   pro=1;
   for m=1:d
   %F0 is a tensor valued function from R^d--R^{d*n_f*d}
       pro=pro*(sum(gw'.*LE(:,NN(i,m)+1).*F0(:,m,l)));
   end
   suma=suma+pro;
end
L(i)=suma;
end

