%For Rastrigin function
function [F0,Nf]=f0parab()
global d gp

F0=ones(length(gp),d,d+1);

for j=1:d
    for k=1:d
        if k==j
        F0(:,k,j)=gp.^2-10.*cos(2.*pi.*gp);
        else
        F0(:,k,j)=1;
        end
    end
end

F0(:,1,d+1)=10*d;
for k=2:d
    F0(:,k,d+1)=1;
end

% separation rank
Nf=d+1;

