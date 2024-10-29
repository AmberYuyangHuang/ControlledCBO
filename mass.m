function [M]=mass(N,NN,d,LELE)

M=sparse(N,N);
for j=1:N
    for i=1:N
        fac=1;
        for k=1:d
            fac=fac*sum(LELE(:,NN(i,k)+1,NN(j,k)+1));
        end
        M(i,j)=fac;
    end
end
