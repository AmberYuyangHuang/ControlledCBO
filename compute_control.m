function [val]=compute_control(x,epsilon,solution)
global NN d Ntotal
LEPP=zeros(Ntotal,d);
for i=1:Ntotal
    for j=1:d
        LEPP(i,j)=legp(x(j),NN(i,j));
        for k=1:d
            if k~=j
                LEPP(i,j)=LEPP(i,j)*leg(x(k),NN(i,k));
            end
        end
    end
end

val=zeros(d,1);
for j=1:d
    for i=1:Ntotal
       val(j)=val(j)+solution(i)*LEPP(i,j);
    end
val(j)=-(1/epsilon)*val(j);
end
