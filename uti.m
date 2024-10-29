function [U1]=uti(sol)

global Ntotal W gamma

U1=zeros(Ntotal,1);

for i=1:Ntotal
    Uk=squeeze(W(:,:,i));
    U1(i)=(0.25/gamma)*(sol'*Uk*sol);
end
