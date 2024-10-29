function [G]=gradit(sol)
global Ntotal W gamma
G=zeros(Ntotal,Ntotal);

for l=1:Ntotal
   for i=1:Ntotal
        vk=squeeze(W(i,:,l));
        G(l,i)=-(0.5/gamma).*vk*sol;
    end
end
        
